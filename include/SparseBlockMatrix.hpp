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
SparseBlockMatrix<BlockType_>::SparseBlockMatrix(const int num_block_rows_,
      const int num_block_cols_, const int block_dim_) {
   _num_block_rows = num_block_rows_;
   _num_block_cols = num_block_cols_;
   _block_dim = block_dim_;
   _total_rows = _num_block_rows * block_dim_;
   _total_cols = _num_block_cols * block_dim_;
   _row_container.resize(num_block_rows_);
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
void SparseBlockMatrix<BlockType_>::resize(const int new_block_rows_, const int new_block_cols_){
   _num_block_rows = new_block_rows_;
   _num_block_cols = new_block_cols_;
   _total_rows = new_block_rows_ *  _block_dim;
   _total_cols = new_block_cols_ * _block_dim;
   _row_container.resize(new_block_rows_);
}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::setBlock(const int r_, const int c_,
      DenseBlock block_){
   if(r_ >= _num_block_rows || c_ >= _num_block_cols){
      cerr << RED << "set Out of bound" << RESET << endl;
      return;
   }

   ColumnsBlockMap& curr_row = _row_container[r_];
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

   if(r_ >= _num_block_rows || c_ >= _num_block_cols){
      cerr << RED << "get Out of bound" << RESET << endl;
      return result;
   }

   const ColumnsBlockMap& curr_row = _row_container[r_];
   typename ColumnsBlockMap::const_iterator it = curr_row.find(c_);
   if(it == curr_row.end())
      return result;
   return it->second;
}

template<typename BlockType_>
bool SparseBlockMatrix<BlockType_>::isNonZeroBlock(const int r_, const int c_) const {
   if(r_ >= _num_block_rows || c_ >= _num_block_cols){
      cerr << RED << "non zero Out of bound" << RESET << endl;
      return false;
   }

   const ColumnsBlockMap& curr_row = _row_container[r_];
   typename ColumnsBlockMap::const_iterator it = curr_row.find(c_);
   if(it == curr_row.end())
      return false;
   return true;
}

template<typename BlockType_>
void SparseBlockMatrix<BlockType_>::printBlock(const int r_, const int c_) const {
   if(r_ >= _num_block_rows || c_ >= _num_block_cols){
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
   block_cholesky_ = SparseBlockMatrix<DenseBlock>(_num_block_rows, _num_block_cols, _block_dim);

   //! Looping over rows
   for (int r = 0; r < _num_block_rows; ++r) {
      const ColumnsBlockMap& curr_row = _row_container[r];
      ColumnsBlockMap& chol_curr_row = block_cholesky_._row_container[r];
      if(curr_row.empty())
         return;
      int starting_col_idx = curr_row.begin()->first;

      //! Looping over cols
      for (int c = starting_col_idx; c <= r; ++c) {
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
void SparseBlockMatrix<BlockType_>::cholesky(SparseBlockMatrix<BlockType_>& block_cholesky_){
   if(_total_cols != _total_rows){
      cerr << BOLDRED << "Error: matrix must be squared" <<
            RESET << endl;
      return;
   }

   block_cholesky_ = SparseBlockMatrix<DenseBlock>(_num_block_rows, _num_block_cols, _block_dim);

   //! Looping over rows
   for (int r = 0; r < _num_block_rows; ++r) {
      const ColumnsBlockMap& curr_row = _row_container[r];
      ColumnsBlockMap& chol_curr_row = block_cholesky_._row_container[r];

      vector<DenseBlock> diagonal_blocks(_num_block_rows);
      vector<DenseBlock> inverse_diagonal_blocks(_num_block_rows);

      if(curr_row.empty())
         return;
      int starting_col_idx = curr_row.begin()->first;

      //! Looping over cols
      for (int c = starting_col_idx; c <= r; ++c) {
         DenseBlock accumulator = DenseBlock::Zero();
         DenseBlock chol_curr_rc_value = DenseBlock::Zero();

         ColumnsBlockMap& chol_upper_row = block_cholesky_._row_container[c];
         accumulator = getBlock(r,c) - scalarProd(chol_curr_row, chol_upper_row, c-1);

         if(r == c){
            //TODO chol_curr_rc_value = sqrt(accumulator);
//            chol_curr_rc_value.setIdentity();
            chol_curr_rc_value = chDecomposeBlock(accumulator);
            diagonal_blocks[r] = chol_curr_rc_value;
            inverse_diagonal_blocks[r] = chol_curr_rc_value.inverse();
         } else {
            chol_curr_rc_value = inverse_diagonal_blocks[c] * accumulator;
         }
         chol_curr_row.insert(make_pair(c, chol_curr_rc_value));
      }
   }
}

template<typename BlockType_>
BlockType_ SparseBlockMatrix<BlockType_>::chDecomposeBlock(const BlockType_& input_){
   DenseBlock lower;
   int dim = input_.rows();
   lower.setZero();
   for (int i = 0; i < dim; i++)
      for (int j = 0; j < (i+1); j++) {
         float s = 0;
         for (int k = 0; k < j; k++)
            s += lower(i, k) * lower(j, k);
         lower(i, j) = (i == j) ?
               sqrt((float)input_(i, j) - s) :
               (1.0 / lower(j,j) * (input_(i,j) - s));
      }
   return lower;
}



template<typename BlockType_>
bool SparseBlockMatrix<BlockType_>::evaluateScalarProdStructure(const ColumnsBlockMap& row_1_,
      const ColumnsBlockMap& row_2_, int max_pos_){
   typename ColumnsBlockMap::const_iterator it_1 = row_1_.begin();
   typename ColumnsBlockMap::const_iterator it_2 = row_2_.begin();

   while(it_1 != row_1_.end() && it_2 != row_2_.end()){
      int col_idx_1 = it_1->first;
      int col_idx_2 = it_2->first;

      if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_)
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

template<typename BlockType_>
BlockType_ SparseBlockMatrix<BlockType_>::scalarProd(const ColumnsBlockMap& row_1_,
      const ColumnsBlockMap& row_2_, int max_pos_){
   typename ColumnsBlockMap::const_iterator it_1 = row_1_.begin();
   typename ColumnsBlockMap::const_iterator it_2 = row_2_.begin();

   DenseBlock result = DenseBlock::Zero();
   while(it_1 != row_1_.end() && it_2 != row_2_.end()){
      int col_idx_1 = it_1->first;
      int col_idx_2 = it_2->first;

      if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_)
         return result;
      if(col_idx_1 == col_idx_2){
         result += it_1->second * it_2->second;
         ++it_1;
         ++it_2;
      }
      else if(col_idx_1 > col_idx_2)
         ++it_2;
      else if(col_idx_1 < col_idx_2)
         ++it_1;
   }
   return result;
}



} /* namespace sparse */
