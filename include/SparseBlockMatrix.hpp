/*
 * SparseBlockMatrix.cpp
 *
 *  Created on: 29/mar/2017
 *      Author: istin
 */
using namespace std;

namespace sparse {


//! ---------------------------------------------------------------- !//
//! --------------------- Methods definitions ---------------------- !//
//! ---------------------------------------------------------------- !//
template<typename BlockType_>
SparseBlockMatrix<BlockType_>::SparseBlockMatrix() {
   _total_rows = 0;
   _total_cols = 0;
   _num_block_rows = 0;
   _num_block_cols = 0;
}

template<typename BlockType_>
SparseBlockMatrix<BlockType_>::SparseBlockMatrix(const int num_rows_,
      const int num_cols_, const int block_dim_) {
   if(num_rows_ % block_dim_ != 0)
      throw std::runtime_error("Error, total dimension and block dimension must be compatible");
   _total_rows = num_rows_;
   _total_cols = num_cols_;
   _block_dim = block_dim_;
   _num_block_rows = _total_rows / block_dim_;
   _num_block_cols = _total_cols / block_dim_;
   _row_container.resize(num_rows_ / _block_dim);
}

template<typename BlockType_>
SparseBlockMatrix<BlockType_>::~SparseBlockMatrix() {

}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::reset(void){
   for(int i = 0; i < _total_rows; ++i){
      _row_container[i].clear();
   }
   _total_rows = 0;
   _total_cols = 0;
   _num_block_rows = 0;
   _num_block_cols = 0;
}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::resize(const int new_rows_, const int new_cols_){
   if(new_rows_ % _block_dim != 0)
      throw std::runtime_error("Error, total dimension and block dimension must be compatible");
   _total_rows = new_rows_;
   _total_cols = new_cols_;
   _num_block_rows = _total_rows / _block_dim;
   _num_block_cols = _total_cols / _block_dim;
   _row_container.resize(new_rows_/_block_dim);
}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::setBlock(const int r_, const int c_,
      DenseBlock block_){
   if (r_%_block_dim != 0 || c_ % _block_dim != 0) {
      cerr << RED << "set bad indexing" << RESET << endl;
      return;
   }
   const int r = r_/_block_dim;
   const int c = c_/_block_dim;
   if(r >= _num_block_rows || c >= _num_block_cols){
      cerr << RED << "set Out of bound" << RESET << endl;
      return;
   }
   ColumnsBlockMap& curr_row = _row_container[r];
   typename ColumnsBlockMap::iterator it = curr_row.find(c_);
   if(it == curr_row.end()){
      if(!block_.isZero()){
         curr_row.insert(std::make_pair(c_, block_));
      }
   } else {
      if(block_.isZero())
         curr_row.erase(it);
      else
         it->second = block_;
   }
   return;
}

template<typename BlockType_>
BlockType_ SparseBlockMatrix<BlockType_>::getBlock(const int r_, const int c_) const {
   DenseBlock result;
   result.setZero();

   if (r_%_block_dim != 0 || c_ % _block_dim != 0) {
      cerr << RED << "get bad indexing" << RESET << endl;
      return result;
   }
   const int r = r_/_block_dim;
   const int c = c_/_block_dim;
   if(r >= _num_block_rows || c >= _num_block_cols){
      cerr << RED << "get Out of bound" << RESET << endl;
      return result;
   }

   const ColumnsBlockMap& curr_row = _row_container[r];
   typename ColumnsBlockMap::const_iterator it = curr_row.find(c_);
   if(it == curr_row.end())
      return result;
   return it->second;
}

template<typename BlockType_>
bool SparseBlockMatrix<BlockType_>::isNonZeroBlock(const int r_, const int c_) const {
   if (r_%_block_dim != 0 || c_ % _block_dim != 0) {
      cerr << RED << "non zero bad indexing" << RESET << endl;
      return false;
   }
   const int r = r_/_block_dim;
   const int c = c_/_block_dim;
   if(r >= _num_block_rows || c >= _num_block_cols){
      cerr << RED << "non zero Out of bound" << RESET << endl;
      return false;
   }

   const ColumnsBlockMap& curr_row = _row_container[r];
   typename ColumnsBlockMap::const_iterator it = curr_row.find(c_);
   if(it == curr_row.end())
      return false;
   return true;
}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::printBlock(const int r_, const int c_) const {
   if (r_%_block_dim != 0 || c_ % _block_dim != 0) {
      return;
   }
   const int r = r_/_block_dim;
   const int c = c_/_block_dim;
   if(r >= _num_block_rows || c >= _num_block_cols){
      cerr << RED << "print Out of bound" << RESET << endl;
      return;
   }

   cout << BOLDWHITE << "Block(" << r_ << "," << c_ << ")" << ":\n" <<
         CYAN << getBlock(r_,c_) << RESET << endl;
}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::evaluateCholeskyStructure(SparseBlockMatrix<BlockType_>& block_cholesky_){
   if(_total_cols != _total_rows){
      cerr << BOLDRED << "Error: matrix must be squared" <<
            RESET << endl;
      return;
   }
   block_cholesky_ = SparseBlockMatrix<DenseBlock>(_total_rows, _total_cols, _block_dim);

   //! Looping over rows
   for (int r = 0; r < _num_block_rows; ++r) {
      const ColumnsBlockMap& curr_row = _row_container[r];
      ColumnsBlockMap& chol_curr_row = block_cholesky_._row_container[r];
      if(curr_row.empty())
         return;
      int starting_col_idx = curr_row.begin()->first;

      //! Looping over cols
      for (int c = starting_col_idx; c < r; ++c) {
         bool not_empty = false;
         if(isNonZeroBlock(r,c)){
            not_empty = true;
         } else {
            ColumnsBlockMap& chol_upper_row = block_cholesky_._row_container[c];
            not_empty = evaluateScalarProdStructure(chol_curr_row, chol_upper_row, c);
         }
         if(not_empty)
            chol_curr_row.insert(make_pair(c, DenseBlock::Ones()));
      }
   }
}

template<typename BlockType_>
bool SparseBlockMatrix<BlockType_>::evaluateScalarProdStructure(const ColumnsBlockMap& row_1_,
      const ColumnsBlockMap& row_2_, int max_pos){
   typename ColumnsBlockMap::const_iterator it_1 = row_1_.begin();
   typename ColumnsBlockMap::const_iterator it_2 = row_2_.begin();

   while(it_1 != row_1_.end() && it_2 != row_2_.end()){
      int col_idx_1 = it_1->first;
      int col_idx_2 = it_2->first;

      if(col_idx_1 > max_pos || col_idx_2 > max_pos)
         return false;
      if(col_idx_1 == col_idx_2)
         return true;
      if(col_idx_1 > col_idx_2)
         ++it_2;
      else if(col_idx_1 < col_idx_2)
         ++it_1;
   }
   return false;
}



} /* namespace sparse */
