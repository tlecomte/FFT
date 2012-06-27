package plugins.praveen.fft;

import cern.colt.function.tdouble.DoubleDoubleFunction;

//function that walks the 2D FFT from JTransforms and fills the sequence data array
//this is the version that does not swap the quadrants
public class DirectAssign2D implements AssignFunction2D {
	public void assign(double[] in, double[] out, int _w, int _h,
			DoubleDoubleFunction function) {
		for (int i = 0; i < in.length/2; i++)
		{
			double real = in[2*i];
			double imag = in[2*i + 1];
			
			out[i] = function.apply(real, imag);
		}	
	}
}
