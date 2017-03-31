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

   sparse::SparseBlockMatrix<Eigen::Matrix2f>* sp_block_mat = new sparse::SparseBlockMatrix<Eigen::Matrix2f>(3,3,2);
   sparse::SparseBlockMatrix<Eigen::Matrix2f> chol;
   Eigen::Matrix2f a;
/*
   1.3137    1.9933         0         0    1.4421    1.7743
   1.9933    3.1610         0         0    2.3542    2.6394
        0         0    0.2354    0.3563         0         0
        0         0    0.3563    0.5396         0         0
   1.4421    2.3542         0         0    2.1761    2.1368
   1.7743    2.6394         0         0    2.1368    2.6851
/**/

   a << 1.3137,   1.9933,
         1.9933, 3.1610;
   sp_block_mat->setBlock(0,0,a);

   a << 1.4421, 1.7743,
         2.3542, 2.6394;
   sp_block_mat->setBlock(0,2,a);

   a << 0.2354, 0.3563,
         0.3563, 0.5396;
   sp_block_mat->setBlock(1,1,a);

   a << 1.4421, 2.3542,
         1.7743, 2.6394;
   sp_block_mat->setBlock(2,0,a);

   a << 2.1761, 2.1368,
         2.1368, 2.6851;
   sp_block_mat->setBlock(2,2,a);

//   sp_block_mat->evaluateCholeskyStructure(chol);
   sp_block_mat->cholesky(chol);




   for (int i = 0; i < sp_block_mat->numRows(); ++i) {
      for (int j = 0; j < sp_block_mat->numCols(); ++j) {
//         sp_block_mat->printBlock(i,j);
         chol.printBlock(i,j);
      }
   }
/**/

   return 0;
}

