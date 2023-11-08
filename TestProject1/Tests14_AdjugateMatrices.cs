﻿using NUnit.Framework;

namespace Maths_Matrices.Tests
{
    [TestFixture]
    public class Tests14_AdjugateMatrices
    {
        [Test]
        public void TestCalculateAdjugateMatrixInstance()
        {
            MatrixFloat m = new MatrixFloat(new[,]
            {
                { 1f, 2f },
                { 3f, 4f }
               
            /* { 1f, 2f,3f },
              { 0f,1f,4f },
              { 5f, 6f, 0f }
            */
            });

            MatrixFloat adjM = m.Adjugate();
            GlobalSettings.DefaultFloatingPointTolerance = 0.001d;
            Assert.AreEqual(new[,]
            {
                { 4f, -3f },
                { -2f, 1f },
                
                /*
                { -24f, 18f,5f },
                { 20f,-15f,-4f },
                { -5f, 4f, 1f}
                */
            }, adjM.ToArray2D());
            GlobalSettings.DefaultFloatingPointTolerance = 0.0d;
            
        }

        [Test]
        public void TestCalculateAdjugateMatrixStatic()
        {
            MatrixFloat m = new MatrixFloat(new[,]
            {
                { 1f, 0f, 5f },
                { 2f, 1f, 6f },
                { 3f, 4f, 0f },
            });

            MatrixFloat adjM = MatrixFloat.Adjugate(m);

            GlobalSettings.DefaultFloatingPointTolerance = 0.001d;
            Assert.AreEqual(new[,]
            {
                { -24f, 18f, 5f },
                { 20f, -15f, -4f },
                { -5f, 4f, 1f },
            }, adjM.ToArray2D());
            GlobalSettings.DefaultFloatingPointTolerance = 0.088d;
        }
    }
}