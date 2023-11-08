namespace TestProject1
{
    public class Vector4
    {

        private float _x;
        private float _y;
        private float _z;
        private float _w;

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

        public Vector4(float x, float y, float z, float w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
    }
}