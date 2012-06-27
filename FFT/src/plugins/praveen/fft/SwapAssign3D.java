package plugins.praveen.fft;

import plugins.praveen.fft.fft.ApplyFunction;

//function that walks the 3D FFT from JTransforms and fills the sequence data array
//this is the version that swaps the quadrants
public class SwapAssign3D implements AssignFunction3D {
	public void assign(double[] in, double[][][] out, int _w, int _h, int _z, int c,
			ApplyFunction function) {
		int wc = (int) Math.ceil(_w/2);
		int hc = (int) Math.ceil(_h/2);
		int zc = (int) Math.ceil(_z/2);
		
		for(int k = 0; k < zc+1; k++)
		{			
			for(int x = 0; x < wc+1; x++)
			{
				for(int y = 0; y < hc+1; y++)
				{
					double real = in[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0];
					double imag = in[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0];
					double imag = in[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
			}
			for(int x = wc+1; x < _w; x++)
			{
				for(int y = 0; y < hc+1; y++)
				{
					double real = in[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0];
					double imag = in[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0];
					double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
			}
		}
		for(int k = zc+1; k < _z; k++)
		{			
			for(int x = 0; x < wc+1; x++)
			{
				for(int y = 0; y < hc+1; y++)
				{
					double real = in[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
					double imag = in[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
					double imag = in[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
			}
			for(int x = wc+1; x < _w; x++)
			{
				for(int y = 0; y < hc+1; y++)
				{
					double real = in[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
					double imag = in[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
					double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
					out[k][c][x + _w*y] = function.apply(real, imag);
				}
			}
		}
	}
}
