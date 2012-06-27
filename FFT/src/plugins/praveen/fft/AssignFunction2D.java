package plugins.praveen.fft;

import cern.colt.function.tdouble.DoubleDoubleFunction;

public interface AssignFunction2D {
	void assign(double[] in, double[] out, int _w, int _h,
			DoubleDoubleFunction Function);
}
