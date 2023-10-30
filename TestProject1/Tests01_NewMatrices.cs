using System;
using System.Collections.Generic;
using System.Diagnostics;
using NUnit.Framework;

namespace Maths_Matrices.Tests
{
    public class MatrixInt
    {
        private int _lines;

        public int Lines
        {
            get => _lines;
            private set => _lines = value;
        }

        private int _columns;

        public int Columns
        {
            get => _columns;
            private set => _columns = value;
        }

        private int[,] matrix;

        public int[,] Matrix
        {
            get => matrix;
            set => matrix = value;
        }

        public int this[int i, int i1]
        {
            get { return Matrix[i, i1]; }
            set { Matrix[i, i1] = value; }
        }

        public MatrixInt(int lines, int columns)
        {
            Matrix = new int[lines, columns];
            Lines = lines;
            Columns = columns;
        }

        public MatrixInt(int[,] matrix)
        {
            Matrix = matrix;
            Lines = matrix.GetLength(0);
            Columns = matrix.GetLength(1);
        }

        public MatrixInt(MatrixInt copy)
        {
            Matrix = new int[copy.Matrix.GetLength(0), copy.Matrix.GetLength(1)];

            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    Matrix[i, j] = copy.Matrix[i, j];
                }
            }

            Lines = copy.Lines;
            Columns = copy.Columns;
        }

        public int[,] ToArray2D()
        {
            return Matrix;
        }


        public static MatrixInt Identity(int size)
        {
            MatrixInt matrixIdentity = new MatrixInt(size, size);

            int identityIndex = 0;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (j == identityIndex)
                    {
                        matrixIdentity.Matrix[i, j] = 1;
                    }
                }

                identityIndex++;
            }

            return matrixIdentity;
        }


        public bool IsIdentity()
        {
            if (Lines != Columns)
            {
                return false;
            }

            for (int i = 0; i < Lines; i++)
            {
                for (int j = 0; j < Lines; j++)
                {
                    if (i == j)
                    {
                        if (Matrix[i, j] != 1)
                        {
                            return false;
                        }
                    }
                    else
                    {
                        if (Matrix[i, j] != 0)
                        {
                            return false;
                        }
                    }
                }
            }

            return true;
        }


        private bool IsIdentityDirty()
        {
            if (Lines != Columns)
            {
                return false;
            }

            int[] identitiesIndex = new int[Matrix.GetLength(0)];
            int indexEmpty = 0;

            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    if (Matrix[i, j] == 1)
                    {
                        if (indexEmpty >= identitiesIndex.Length)
                        {
                            return false;
                        }

                        identitiesIndex[indexEmpty] = j;
                        indexEmpty++;
                    }
                    else if (Matrix[i, j] != 0)
                    {
                        return false;
                    }
                }
            }

            for (int i = 0; i < identitiesIndex.Length; i++)
            {
                if (i > 0)
                {
                    if (identitiesIndex[i] - identitiesIndex[i - 1] != 1)
                    {
                        return false;
                    }
                }
            }

            return true;
        }


        public void Multiply(int factor)
        {
            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    Matrix[i, j] *= factor;
                }
            }
        }

        public MatrixInt Multiply(MatrixInt factor)
        {
            if (Columns != factor.Lines)
            {
                throw new MatrixSumException();
            }


            /*for (int i = 0; i < factor.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < factor.Matrix.GetLength(1); j++)
                {
                    Matrix[i,j] = Matrix[i,j] * factor.Matrix[]
                }
            }
            */


            return this;
        }


        public static MatrixInt Multiply(MatrixInt matrixClass, int factor)
        {
            MatrixInt newMatrix = new MatrixInt(matrixClass);

            for (int i = 0; i < newMatrix.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < newMatrix.Matrix.GetLength(1); j++)
                {
                    newMatrix.Matrix[i, j] *= factor;
                }
            }

            return newMatrix;
        }

        public void Add(MatrixInt a)
        {
            if (a.Columns != Columns || a.Lines != Lines)
            {
                throw new MatrixSumException();
            }

            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    Matrix[i, j] += a.Matrix[i, j];
                }
            }
        }

        public static MatrixInt Add(MatrixInt a, MatrixInt b)
        {
            if (a.Columns != b.Columns || a.Lines != b.Lines)
            {
                throw new MatrixSumException();
            }

            MatrixInt newMatrix = new MatrixInt(a);

            for (int i = 0; i < a.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < a.Matrix.GetLength(1); j++)
                {
                    newMatrix.Matrix[i, j] = a.Matrix[i, j] + b.Matrix[i, j];
                }
            }

            return newMatrix;
        }

        public static MatrixInt operator +(MatrixInt a, MatrixInt b) => Add(a, b);
        public static MatrixInt operator *(MatrixInt a, int b) => Multiply(a, b);
        public static MatrixInt operator *(int b, MatrixInt a) => Multiply(a, b);

        public static MatrixInt operator -(MatrixInt a, MatrixInt b)
        {
            if (a.Columns != b.Columns || a.Lines != b.Lines)
            {
                throw new MatrixSumException();
            }

            MatrixInt newMatrix = new MatrixInt(a);

            for (int i = 0; i < a.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < a.Matrix.GetLength(1); j++)
                {
                    newMatrix[i, j] = a.Matrix[i, j] - b.Matrix[i, j];
                }
            }

            return newMatrix;
        }


        public static MatrixInt operator -(MatrixInt a)
        {
            MatrixInt newMatrix = new MatrixInt(a);

            for (int i = 0; i < newMatrix.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < newMatrix.Matrix.GetLength(1); j++)
                {
                    newMatrix.Matrix[i, j] = -a.Matrix[i, j];
                }
            }

            return newMatrix;
        }

        public MatrixInt Transpose()
        {
            throw new NotImplementedException();
        }

        /*  public static MatrixInt Transpose(MatrixInt a)
          {
              return new MatrixInt();
          }
          */
        public static MatrixInt GenerateAugmentedMatrix(MatrixInt mTransfo, MatrixInt mCol)
        {
            MatrixInt augmentedMatrix = new MatrixInt(mTransfo.Lines, mTransfo.Columns + 1);

            int index = 0;
            for (int i = 0; i < augmentedMatrix.Lines; i++)
            {
                for (int j = 0; j < augmentedMatrix.Columns; j++)
                {
                    if (j < augmentedMatrix.Columns - 1)
                    {
                        augmentedMatrix[i, j] = mTransfo[i, j];
                    }
                    else
                    {
                        augmentedMatrix[i, j] =
                            mCol[index, 0]; // Peut avoir plusieurs colonnes cas a gérer dans les prochains test
                        index++;
                    }
                }
            }


            return augmentedMatrix;
        }

        public (MatrixInt a, MatrixInt b) Split(int split)
        {
            MatrixInt firstMatrix = new MatrixInt(Lines, (split + 1));
            MatrixInt secondMatrix = new MatrixInt(Lines, Columns - (split + 1));

            for (int i = 0; i < firstMatrix.Lines; i++)
            {
                for (int j = 0; j < firstMatrix.Columns; j++)
                {
                    firstMatrix[i, j] = Matrix[i, j];
                }

                for (int j = 0; j < secondMatrix.Columns; j++)
                {
                    secondMatrix[i, j] = Matrix[i, j + firstMatrix.Columns];
                }
            }


            return (firstMatrix, secondMatrix);
        }
    }

    public class MatrixFloat
    {
        private int _lines;

        public int Lines
        {
            get => _lines;
            private set => _lines = value;
        }

        private int _columns;

        public int Columns
        {
            get => _columns;
            private set => _columns = value;
        }

        private float[,] matrix;

        public float[,] Matrix
        {
            get => matrix;
            set => matrix = value;
        }

        public float this[int i, int i1]
        {
            get { return Matrix[i, i1]; }
            set { Matrix[i, i1] = value; }
        }

        public MatrixFloat(int lines, int columns)
        {
            Matrix = new float[lines, columns];
            Lines = lines;
            Columns = columns;
        }

        public MatrixFloat(float[,] matrix)
        {
            Matrix = matrix;
            Lines = matrix.GetLength(0);
            Columns = matrix.GetLength(1);
        }

        public MatrixFloat(MatrixFloat copy)
        {
            Matrix = new float[copy.Matrix.GetLength(0), copy.Matrix.GetLength(1)];

            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    Matrix[i, j] = copy.Matrix[i, j];
                }
            }

            Lines = copy.Lines;
            Columns = copy.Columns;
        }

        public float[,] ToArray2D()
        {
            return Matrix;
        }


        public static MatrixFloat Identity(int size)
        {
            MatrixFloat matrixIdentity = new MatrixFloat(size, size);

            int identityIndex = 0;
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (j == identityIndex)
                    {
                        matrixIdentity.Matrix[i, j] = 1;
                    }
                }

                identityIndex++;
            }

            return matrixIdentity;
        }


        public bool IsIdentity()
        {
            if (Lines != Columns)
            {
                return false;
            }

            for (int i = 0; i < Lines; i++)
            {
                for (int j = 0; j < Lines; j++)
                {
                    if (i == j)
                    {
                        if (Matrix[i, j] != 1)
                        {
                            return false;
                        }
                    }
                    else
                    {
                        if (Matrix[i, j] != 0)
                        {
                            return false;
                        }
                    }
                }
            }

            return true;
        }


        public void Multiply(int factor)
        {
            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    Matrix[i, j] *= factor;
                }
            }
        }

        public MatrixFloat Multiply(MatrixFloat factor)
        {
            if (Columns != factor.Lines)
            {
                throw new MatrixSumException();
            }


            /*for (int i = 0; i < factor.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < factor.Matrix.GetLength(1); j++)
                {
                    Matrix[i,j] = Matrix[i,j] * factor.Matrix[]
                }
            }
            */


            return this;
        }


        public static MatrixFloat Multiply(MatrixFloat matrixClass, int factor)
        {
            MatrixFloat newMatrix = new MatrixFloat(matrixClass);

            for (int i = 0; i < newMatrix.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < newMatrix.Matrix.GetLength(1); j++)
                {
                    newMatrix.Matrix[i, j] *= factor;
                }
            }

            return newMatrix;
        }

        public void Add(MatrixFloat a)
        {
            if (a.Columns != Columns || a.Lines != Lines)
            {
                throw new MatrixSumException();
            }

            for (int i = 0; i < Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < Matrix.GetLength(1); j++)
                {
                    Matrix[i, j] += a.Matrix[i, j];
                }
            }
        }

        public static MatrixFloat Add(MatrixFloat a, MatrixFloat b)
        {
            if (a.Columns != b.Columns || a.Lines != b.Lines)
            {
                throw new MatrixSumException();
            }

            MatrixFloat newMatrix = new MatrixFloat(a);

            for (int i = 0; i < a.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < a.Matrix.GetLength(1); j++)
                {
                    newMatrix.Matrix[i, j] = a.Matrix[i, j] + b.Matrix[i, j];
                }
            }

            return newMatrix;
        }

        public static MatrixFloat operator +(MatrixFloat a, MatrixFloat b) => Add(a, b);
        public static MatrixFloat operator *(MatrixFloat a, int b) => Multiply(a, b);
        public static MatrixFloat operator *(int b, MatrixFloat a) => Multiply(a, b);

        public static MatrixFloat operator -(MatrixFloat a, MatrixFloat b)
        {
            if (a.Columns != b.Columns || a.Lines != b.Lines)
            {
                throw new MatrixSumException();
            }

            MatrixFloat newMatrix = new MatrixFloat(a);

            for (int i = 0; i < a.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < a.Matrix.GetLength(1); j++)
                {
                    newMatrix[i, j] = a.Matrix[i, j] - b.Matrix[i, j];
                }
            }

            return newMatrix;
        }


        public static MatrixFloat operator -(MatrixFloat a)
        {
            MatrixFloat newMatrix = new MatrixFloat(a);

            for (int i = 0; i < newMatrix.Matrix.GetLength(0); i++)
            {
                for (int j = 0; j < newMatrix.Matrix.GetLength(1); j++)
                {
                    newMatrix.Matrix[i, j] = -a.Matrix[i, j];
                }
            }

            return newMatrix;
        }

        public MatrixFloat Transpose()
        {
            throw new NotImplementedException();
        }

        /*  public static MatrixInt Transpose(MatrixInt a)
          {
              return new MatrixInt();
          }
          */
        public static MatrixFloat GenerateAugmentedMatrix(MatrixFloat mTransfo, MatrixFloat mCol)
        {
            MatrixFloat augmentedMatrix = new MatrixFloat(mTransfo.Lines, mTransfo.Columns + 1);

            int index = 0;
            for (int i = 0; i < augmentedMatrix.Lines; i++)
            {
                for (int j = 0; j < augmentedMatrix.Columns; j++)
                {
                    if (j < augmentedMatrix.Columns - 1)
                    {
                        augmentedMatrix[i, j] = mTransfo[i, j];
                    }
                    else
                    {
                        augmentedMatrix[i, j] =
                            mCol[index, 0]; // Peut avoir plusieurs colonnes cas a gérer dans les prochains test
                        index++;
                    }
                }
            }


            return augmentedMatrix;
        }

        public (MatrixFloat a, MatrixFloat b) Split(int split)
        {
            MatrixFloat firstMatrix = new MatrixFloat(Lines, (split + 1));
            MatrixFloat secondMatrix = new MatrixFloat(Lines, Columns - (split + 1));

            for (int i = 0; i < firstMatrix.Lines; i++)
            {
                for (int j = 0; j < firstMatrix.Columns; j++)
                {
                    firstMatrix[i, j] = Matrix[i, j];
                }

                for (int j = 0; j < secondMatrix.Columns; j++)
                {
                    secondMatrix[i, j] = Matrix[i, j + firstMatrix.Columns];
                }
            }


            return (firstMatrix, secondMatrix);
        }
    }

    public class MatrixRowReductionAlgorithm
    {
        public static (MatrixFloat m1, MatrixFloat m2) Apply(MatrixFloat m1, MatrixFloat m2)
        {
            MatrixFloat reducMatrix = MatrixFloat.GenerateAugmentedMatrix(m1, m2);


            for (int i = 0; i < m1.Lines; i++)
            {
                for (int j = i; j < m1.Columns; j++)
                {
                    (int lineIndex, float k) = FindHigherK(reducMatrix, i, j);

                    if (AllIsZero(reducMatrix, i))
                    {
                        continue;
                    }

                    if (lineIndex != i)
                    {
                        MatrixElementaryOperations.SwapLines(reducMatrix, lineIndex, i);
                    }

                    MatrixElementaryOperations.MultiplyLine(reducMatrix, i, 1 / reducMatrix[i, j]);

                    for (int l = 0; l < reducMatrix.Lines; l++)
                    {
                        if (l != i)
                        {
                            float targetValueLine = reducMatrix[l, i];
                            
                            for (int r = 0; r < reducMatrix.Columns; r++)
                            {
                                float factor = reducMatrix[j, r];
                                reducMatrix[l, r] -= targetValueLine * factor;
                            }
                        }
                    }

                    break;
                }
            }
            
            return reducMatrix.Split(2);
        }

        private static (int, float) FindHigherK(MatrixFloat matrixRef, int line, int column)
        {
            float highestValue = matrixRef[line, column];
            int lineIndex = 0;

            for (int k = 0; k < matrixRef.Lines; k++)
            {
                if (k >= line)
                {
                    if (matrixRef[k, column] >= highestValue)
                    {
                        highestValue = matrixRef[k, column];
                        lineIndex = k;
                    }
                }
            }

            return (lineIndex, highestValue);
        }

        private static bool AllIsZero(MatrixFloat matrixRef, int line)
        {
            for (int i = 0; i < matrixRef.Columns; i++)
            {
                if (matrixRef[line, i] != 0)
                {
                    return false;
                }
            }

            return true;
        }
    }

    public class MatrixElementaryOperations
    {
        public static void SwapLines(MatrixInt matrixRef, int line1, int line2)
        {
            MatrixInt newMatrix = new MatrixInt(matrixRef);

            for (int i = 0; i < matrixRef.Columns; i++)
            {
                matrixRef[line1, i] = newMatrix[line2, i];
                matrixRef[line2, i] = newMatrix[line1, i];
            }
        }

        public static void SwapLines(MatrixFloat matrixRef, int line1, int line2)
        {
            MatrixFloat newMatrix = new MatrixFloat(matrixRef);

            for (int i = 0; i < matrixRef.Columns; i++)
            {
                matrixRef[line1, i] = newMatrix[line2, i];
                matrixRef[line2, i] = newMatrix[line1, i];
            }
        }


        public static void SwapColumns(MatrixInt matrixRef, int column1, int column2)
        {
            MatrixInt newMatrix = new MatrixInt(matrixRef);

            for (int i = 0; i < matrixRef.Lines; i++)
            {
                matrixRef[i, column1] = newMatrix[i, column2];
                matrixRef[i, column2] = newMatrix[i, column1];
            }
        }

        public static void MultiplyLine(MatrixFloat matrixRef, int line, float factor)
        {
            if (factor == 0)
            {
                throw new MatrixScalarZeroException();
            }

            for (int i = 0; i < matrixRef.Columns; i++)
            {
                matrixRef[line, i] *= factor;
            }
        }

        public static void MultiplyLine(MatrixInt matrixRef, int line, int factor)
        {
            if (factor == 0)
            {
                throw new MatrixScalarZeroException();
            }

            for (int i = 0; i < matrixRef.Columns; i++)
            {
                matrixRef[line, i] *= factor;
            }
        }

        public static void MultiplyColumn(MatrixInt matrixRef, int column, int factor)
        {
            if (factor == 0)
            {
                throw new MatrixScalarZeroException();
            }

            for (int i = 0; i < matrixRef.Lines; i++)
            {
                matrixRef[i, column] *= factor;
            }
        }

        public static void AddLineToAnother(MatrixInt matrixRef, int lineToAdd, int lineResult, int factor)
        {
            for (int i = 0; i < matrixRef.Columns; i++)
            {
                matrixRef[lineResult, i] += matrixRef[lineToAdd, i] * factor;
            }
        }

        public static void AddLineToAnother(MatrixFloat matrixRef, int lineToAdd, int lineResult, float factor)
        {
            for (int i = 0; i < matrixRef.Columns; i++)
            {
                matrixRef[lineResult, i] -= matrixRef[lineToAdd, i] * factor;
            }
        }

        public static void AddColumnToAnother(MatrixInt matrixRef, int columnToAdd, int columnResult, int factor)
        {
            for (int i = 0; i < matrixRef.Lines; i++)
            {
                matrixRef[i, columnResult] = matrixRef[i, columnToAdd] * factor + matrixRef[i, columnResult];
            }
        }
    }

    [TestFixture]
    public class Tests01_NewMatrices
    {
        [Test]
        public void TestNewEmptyMatrices()
        {
            MatrixInt m1 = new MatrixInt(3, 2);
            Assert.AreEqual(new[,]
            {
                { 0, 0 },
                { 0, 0 },
                { 0, 0 }
            }, m1.ToArray2D());
            Assert.AreEqual(3, m1.Lines);
            Assert.AreEqual(2, m1.Columns);

            MatrixInt m2 = new MatrixInt(2, 3);
            Assert.AreEqual(new[,]
            {
                { 0, 0, 0 },
                { 0, 0, 0 },
            }, m2.ToArray2D());
            Assert.AreEqual(2, m2.Lines);
            Assert.AreEqual(3, m2.Columns);
        }

        [Test]
        public void TestNewMatricesFrom2DArray()
        {
            //See GetLength documentation to retrieve length of a multi-dimensional array
            //https://docs.microsoft.com/en-us/dotnet/api/system.array.getlength
            MatrixInt m = new MatrixInt(new[,]
                {
                    { 1, 2, 3 },
                    { 4, 5, 6 },
                    { 7, 8, 9 },
                }
            );
            Assert.AreEqual(3, m.Lines);
            Assert.AreEqual(3, m.Columns);

            //See Indexers Documentation =>
            //https://docs.microsoft.com/fr-fr/dotnet/csharp/programming-guide/indexers/
            Assert.AreEqual(1, m[0, 0]);
            Assert.AreEqual(2, m[0, 1]);
            Assert.AreEqual(3, m[0, 2]);
            Assert.AreEqual(4, m[1, 0]);
            Assert.AreEqual(5, m[1, 1]);
            Assert.AreEqual(6, m[1, 2]);
            Assert.AreEqual(7, m[2, 0]);
            Assert.AreEqual(8, m[2, 1]);
            Assert.AreEqual(9, m[2, 2]);
        }
    }
}