#include <iostream>
#include <Eigen/Core>

#include "SparseMatrix.h"
#include "SparseBlockMatrix.h"

int main(int argc, char const *argv[])
{
   /*
   sparse::SparseMatrix* sp_mat = new sparse::SparseMatrix();
   sp_mat->loadFromTXT("../data/sparse_data.txt");
//   sp_mat->loadFromTXT("../data/small_sparse_matrix.txt");
   sparse::SparseMatrix chol(0,0);
   sp_mat->cholesky(chol);
/**/

   Eigen::Matrix3f m;
   m.setZero();
   //   m->setZero();

   /*
   sparse::SparseBlockMatrix<Eigen::Matrix3f>* sp_block_mat = new sparse::SparseBlockMatrix<Eigen::Matrix3f>(12,12,3);
   sp_block_mat->setBlock(0,0,m);

   m.setIdentity();
   sp_block_mat->setBlock(3,0,m);

   m.setRandom();
   sp_block_mat->setBlock(6,0,Eigen::Matrix3f::Ones());

   sparse::SparseBlockMatrix<Eigen::Matrix3f> chol;
   sp_block_mat->evaluateCholeskyStructure(chol);
/**/

   sparse::SparseBlockMatrix<Eigen::Matrix2f>* sp_block_mat = new sparse::SparseBlockMatrix<Eigen::Matrix2f>(6,6,2);
   sparse::SparseBlockMatrix<Eigen::Matrix2f> chol;
   Eigen::Matrix2f a;

   a << 2.744007,  1.919497,
         1.919497, 3.191362;
   sp_block_mat->setBlock(0,0,a);

   a << 1.806722,  2.258476,
         2.390472, 2.448052;
   sp_block_mat->setBlock(0,2,a);

   a << 2.123150,  0.823146,
         1.962895, 1.338299;
   sp_block_mat->setBlock(0,4,a);

   a << 1.806722,  2.390472,
         2.258476, 2.448052;
   sp_block_mat->setBlock(2,0,a);

   a << 2.828793, 2.967253,
         2.967253, 3.572853;
   sp_block_mat->setBlock(2,2,a);

   a << 2.101034, 1.714458,
         2.482527, 1.624698;
   sp_block_mat->setBlock(2,4,a);

   a << 2.123150,  1.962895,
         0.823146, 1.338299;
   sp_block_mat->setBlock(4,0,a);

   a << 2.101034,  2.482527,
         1.714458, 1.624698;
   sp_block_mat->setBlock(4,2,a);

   a << 2.199896,  0.931872,
         0.931872, 1.265850;
   sp_block_mat->setBlock(4,4,a);

   sp_block_mat->evaluateCholeskyStructure(chol);

   for (int i = 0; i < sp_block_mat->numRows(); ++i) {
      for (int j = 0; j < sp_block_mat->numCols(); ++j) {
//         sp_block_mat->printBlock(i,j);
         chol.printBlock(i,j);
      }
   }
   /**/

   return 0;
}

