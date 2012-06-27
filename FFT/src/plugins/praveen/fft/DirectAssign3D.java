package plugins.praveen.fft;

import plugins.praveen.fft.fft.ApplyFunction;

//function that walks the 3D FFT from JTransforms and fills the sequence data array
//this is the version that does not swap the quadrants
public class DirectAssign3D implements AssignFunction3D {
	public void assign(double[] in, double[][][] out, int _w, int _h, int _z, int c,
			ApplyFunction function) {
		for(int k = 0; k < _z; k++)
		{			
			for(int x = 0; x < _w; x++)
			{
				for(int y = 0; y < _h; y++)
				{
					double real = in[(x + (y * _w) + (k * _w * _h))*2 + 0];
					double imag = in[(x + (y * _w) + (k * _w * _h))*2 + 1];					
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
			}
		}
	}	
}
