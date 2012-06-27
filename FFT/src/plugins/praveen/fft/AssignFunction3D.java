package plugins.praveen.fft;

import plugins.praveen.fft.fft.ApplyFunction;

public interface AssignFunction3D {
	void assign(double[] in, double[][][] out, int _w, int _h, int _z, int c,
			ApplyFunction function);
}
