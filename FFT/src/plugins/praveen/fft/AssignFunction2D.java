package plugins.praveen.fft;

import plugins.praveen.fft.fft.ApplyFunction;

public interface AssignFunction2D {
	void assign(double[] in, double[] out, int _w, int _h,
			ApplyFunction Function);
}
