using Maths_Matrices.Tests;

namespace TestProject1
{
    public class Vector4
    {

        private float _x;
        private float _y;
        private float _z;
        private float _w;
        
        private MatrixFloat _matrix;

        public float x
        {
            get => _x;
            set => _x = value;
        }

        public float y
        {
            get => _y;
            set => _y = value;
        }

        public float z
        {
            get => _z;
            set => _z = value;
        }

        public float w
        {
            get => _w;
            set => _w = value;
        }

        public MatrixFloat Matrix
        {
            get => _matrix;
            set => _matrix = value;
        }

        public Vector4(float x, float y, float z, float w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;

            _matrix = new MatrixFloat(4,1);
            _matrix[0, 0] = x;
            _matrix[1, 0] = y;
            _matrix[2, 0] = z;
            _matrix[3, 0] = w;
        }
        
        public static  Vector4 operator *(Vector4 vector,MatrixFloat matrix) => MatrixFloat.Multiply(matrix,vector.Matrix).ToVector4();
        public static  Vector4 operator *(MatrixFloat matrix,Vector4 vector) => MatrixFloat.Multiply(matrix,vector.Matrix).ToVector4();
        
        }
}