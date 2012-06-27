package plugins.praveen.fft;

import cern.colt.function.tdouble.DoubleDoubleFunction;

// function that walks the 2D FFT from JTransforms and fills the sequence data array
// this is the version that swaps the quadrants
public class SwapAssign2D implements AssignFunction2D{
	public void assign(double[] in, double[] out, int _w, int _h,
			DoubleDoubleFunction function)
	{
		int wc = (int) Math.ceil(_w/2);
		int hc = (int) Math.ceil(_h/2);
		
		for(int x = 0; x < (wc+1); x++)
		{
			for(int y = 0; y < (hc+1); y++)
			{
				double real = in[((wc-x) + (hc-y) * _w)*2 + 0];
				double imag = in[((wc-x) + (hc-y) * _w)*2 + 1];
				
				out[x + _w*y] = function.apply(real, imag);
			}
			for(int y = hc+1; y < _h; y++)
			{
				double real = in[((wc-x) + (_h+(hc-y)) * _w)*2 + 0];
				double imag = in[((wc-x) + (_h+(hc-y)) * _w)*2 + 1];
				
				out[x + _w*y] = function.apply(real, imag);
			}
	
		}
		for(int x = (wc+1); x < _w; x++)
		{
			for(int y = 0; y < (hc+1); y++)
			{
				double real = in[((_w+(wc-x)) + (hc-y) * _w)*2 + 0];
				double imag = in[((_w+(wc-x)) + (hc-y) * _w)*2 + 1];
				
				out[x + _w*y] = function.apply(real, imag);
			}
			for(int y = hc+1; y < _h; y++)
			{
				double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w)*2 + 0];
				double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w)*2 + 1];
				
				out[x + _w*y] = function.apply(real, imag);
			}
		}
	}
}
