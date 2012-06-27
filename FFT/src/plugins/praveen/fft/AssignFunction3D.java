package plugins.praveen.fft;

import cern.colt.function.tdouble.DoubleDoubleFunction;

public interface AssignFunction3D {
	void assign(double[] in, double[][][] out, int _w, int _h, int _z, int c,
			DoubleDoubleFunction function);
}
