package plugins.praveen.fft;

import cern.colt.function.tdouble.DoubleDoubleFunction;

public class ComplexFunctions {
	public static class Real implements DoubleDoubleFunction
	{
		public double apply(double real, double imag)
		{
			return real;
		}
	};
	
	public static class Imag implements DoubleDoubleFunction
	{
		public double apply(double real, double imag)
		{
			return imag;
		}
	};
	
	public static class Magnitude implements DoubleDoubleFunction
	{
		public double apply(double real, double imag)
		{
			return Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
		}
	};

	public static class Angle implements DoubleDoubleFunction
	{
		public double apply(double real, double imag)
		{
			return Math.atan2(imag, real);
		}
	};
}
