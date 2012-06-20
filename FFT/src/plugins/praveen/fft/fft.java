package plugins.praveen.fft;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.vars.lang.VarSequence;

public class fft extends EzPlug implements Block {

	enum FFTDims {
		FFT_2D("2D"), FFT_3D("3D");
		private String stringValue;
		FFTDims(String s) { stringValue = s; }
		public String toString() { return stringValue; }
	}
	
	enum FFTOutputType {
		MAGNITUDE_PHASE("Magnitude/Phase Pair"), REAL_IMAG("Real/Imaginary Pair");
		private String stringValue;
		FFTOutputType(String s) { stringValue = s; }
		public String toString() { return stringValue; }
	}
	
	EzVarSequence input = new EzVarSequence("Input");
	EzVarEnum<FFTDims> ndims =  new EzVarEnum<FFTDims>("Type", FFTDims.values(), 0);
	EzVarEnum<FFTOutputType> outputType =  new EzVarEnum<FFTOutputType>("Output as", FFTOutputType.values(), 0);
	EzVarBoolean	swap = new EzVarBoolean("Swap Quadrants?", false);
	
	VarSequence fSequenceVar = new VarSequence("FFT sequence", null);

	@Override
	protected void initialize() {
		super.addEzComponent(input);
		super.addEzComponent(ndims);
		super.addEzComponent(outputType);
		super.addEzComponent(swap);
		super.setTimeDisplay(true);
	}
	
	// declare ourself to Blocks
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(input.getVariable());
		inputMap.add(ndims.getVariable());
		inputMap.add(outputType.getVariable());
		inputMap.add(swap.getVariable());
	}

	// declare ourself to Blocks
	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add(fSequenceVar);
	}

	@Override
	protected void execute() {
		Sequence sequence = input.getValue();
		Sequence fSequence = null;

		if(ndims.getValue()==FFTDims.FFT_2D)		
			fSequence = FFT_2D(sequence, swap.getValue(), outputType.getValue());	
		else
			fSequence = FFT_3D(sequence, swap.getValue(), outputType.getValue());
		
		if (!isHeadLess()) {
			addSequence(fSequence);
		}
		
		fSequenceVar.setValue(fSequence);
	}
	
	interface ApplyFunction {
		double apply(double real, double imag);
	}
	
	ApplyFunction realApplyFunction = new ApplyFunction()
	{
		public double apply(double real, double imag)
		{
			return real;
		}
	};
	
	ApplyFunction imagApplyFunction = new ApplyFunction()
	{
		public double apply(double real, double imag)
		{
			return imag;
		}
	};
	
	ApplyFunction magnitudeApplyFunction = new ApplyFunction()
	{
		public double apply(double real, double imag)
		{
			return Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
		}
	};

	ApplyFunction angleApplyFunction = new ApplyFunction()
	{
		public double apply(double real, double imag)
		{
			return Math.atan2(imag, real);
		}
	};

	
	interface AssignFunction3D {
		void assign(double[] in, double[][][] out, int _w, int _h, int _z,
				ApplyFunction channel0Function,
				ApplyFunction channel1Function);
	}
	
	// function that walks the 3D FFT from JTransforms and fills the sequence data array
	// this is the version that does not swap the quadrants
	AssignFunction3D directAssign3D = new AssignFunction3D()
	{
		public void assign(double[] in, double[][][] out, int _w, int _h, int _z,
				ApplyFunction channel0Function,
				ApplyFunction channel1Function) {
			for(int k = 0; k < _z; k++)
			{			
				for(int x = 0; x < _w; x++)
				{
					for(int y = 0; y < _h; y++)
					{
						double real = in[(x + (y * _w) + (k * _w * _h))*2 + 0];
						double imag = in[(x + (y * _w) + (k * _w * _h))*2 + 1];					
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
				}
			}
		}	
	};

	// function that walks the 3D FFT from JTransforms and fills the sequence data array
	// this is the version that swaps the quadrants
	AssignFunction3D swapAssign3D = new AssignFunction3D()
	{
		public void assign(double[] in, double[][][] out, int _w, int _h, int _z,
				ApplyFunction channel0Function,
				ApplyFunction channel1Function) {
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
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0];
						double imag = in[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
				}
				for(int x = wc+1; x < _w; x++)
				{
					for(int y = 0; y < hc+1; y++)
					{
						double real = in[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
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
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
						double imag = in[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
				}
				for(int x = wc+1; x < _w; x++)
				{
					for(int y = 0; y < hc+1; y++)
					{
						double real = in[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = channel0Function.apply(real, imag);
						out[k][1][x + _w*y] = channel1Function.apply(real, imag);
					}
				}
			}
		}
	};

	private Sequence FFT_3D(Sequence sequence, boolean swap, FFTOutputType outputType) {
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();

		final DoubleFFT_3D fft = new DoubleFFT_3D(_z, _h, _w);
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 3D");

		// allocate the output sequence
		for(int k = 0; k < _z; k++)
		{	
			IcyBufferedImage fImage = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);
			fSequence.setImage(0, k, fImage);			
		}
		
		double[] fArray = new double[_w*_h*_z*2];
		// copy the data in fArray, with proper structure
		for(int k = 0; k < _z; k++)
		{
			Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), 0, fArray, k*_w*_h, _w*_h, sequence.isSignedDataType());
		}

		fft.realForwardFull(fArray);
		
		// direct reference to 3D byte array data [Z][C][XY] for specified t
		double[][][] resultData = fSequence.getDataXYCZAsDouble(0);

		ApplyFunction channel0ApplyFunction = null;
		ApplyFunction channel1ApplyFunction = null;
		if(outputType == FFTOutputType.MAGNITUDE_PHASE)
		{
			channel0ApplyFunction = magnitudeApplyFunction;
			channel1ApplyFunction = angleApplyFunction;
			fSequence.setChannelName(0, "Magnitude");
			fSequence.setChannelName(1, "Phase");
		}
		else
		{
			channel0ApplyFunction = realApplyFunction;
			channel1ApplyFunction = imagApplyFunction;
			fSequence.setChannelName(0, "Real");
			fSequence.setChannelName(1, "Imaginary");
		}
		
		AssignFunction3D assignFunction = null;
		if(!swap)  // No Quadrant swapping. Leave as it is.
		{
			assignFunction = directAssign3D;
		}
		else
		{
			assignFunction = swapAssign3D; // Swap Quadrants
		}
		
		assignFunction.assign(fArray, resultData, _w, _h, _z, channel0ApplyFunction, channel1ApplyFunction);

		fSequence.dataChanged();
		
		return fSequence;
	}
	
	interface AssignFunction2D {
		void assign(double[] in, double[][] out, int _w, int _h,
				ApplyFunction channel0Function,
				ApplyFunction channel1Function);
	}
	
	// function that walks the 2D FFT from JTransforms and fills the sequence data array
	// this is the version that does not swap the quadrants
	AssignFunction2D directAssign = new AssignFunction2D()
	{
		public void assign(double[] in, double[][] out, int _w, int _h,
				ApplyFunction channel0Function,
				ApplyFunction channel1Function) {
			for (int i = 0; i < in.length/2; i++)
			{
				double real = in[2*i];
				double imag = in[2*i + 1];
				
				out[0][i] = channel0Function.apply(real, imag);
				out[1][i] = channel1Function.apply(real, imag);
			}	
		}	
	};
	
	// function that walks the 2D FFT from JTransforms and fills the sequence data array
	// this is the version that swaps the quadrants
	AssignFunction2D swapAssign = new AssignFunction2D() {
		public void assign(double[] in, double[][] out, int _w, int _h,
				ApplyFunction channel0Function,
				ApplyFunction channel1Function)
		{
			int wc = (int) Math.ceil(_w/2);
			int hc = (int) Math.ceil(_h/2);
			
			for(int x = 0; x < (wc+1); x++)
			{
				for(int y = 0; y < (hc+1); y++)
				{
					double real = in[((wc-x) + (hc-y) * _w)*2 + 0];
					double imag = in[((wc-x) + (hc-y) * _w)*2 + 1];
					
					out[0][x + _w*y] = channel0Function.apply(real, imag);
					out[1][x + _w*y] = channel1Function.apply(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((wc-x) + (_h+(hc-y)) * _w)*2 + 0];
					double imag = in[((wc-x) + (_h+(hc-y)) * _w)*2 + 1];
					
					out[0][x + _w*y] = channel0Function.apply(real, imag);
					out[1][x + _w*y] = channel1Function.apply(real, imag);
				}
		
			}
			for(int x = (wc+1); x < _w; x++)
			{
				for(int y = 0; y < (hc+1); y++)
				{
					double real = in[((_w+(wc-x)) + (hc-y) * _w)*2 + 0];
					double imag = in[((_w+(wc-x)) + (hc-y) * _w)*2 + 1];
					
					out[0][x + _w*y] = channel0Function.apply(real, imag);
					out[1][x + _w*y] = channel1Function.apply(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w)*2 + 0];
					double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w)*2 + 1];
					
					out[0][x + _w*y] = channel0Function.apply(real, imag);
					out[1][x + _w*y] = channel1Function.apply(real, imag);
				}
			}
		}
	};

	private Sequence FFT_2D(Sequence sequence, boolean swap, FFTOutputType outputType) 
	{
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 2D");
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();

		final DoubleFFT_2D fft = new DoubleFFT_2D(_h, _w);
		
		ApplyFunction channel0Function = null;
		ApplyFunction channel1Function = null;
		if(outputType == FFTOutputType.MAGNITUDE_PHASE)
		{
			channel0Function = magnitudeApplyFunction;
			channel1Function = angleApplyFunction;
			fSequence.setChannelName(0, "Magnitude");
			fSequence.setChannelName(1, "Phase");
		}
		else // Real/Imaginary Pair
		{
			channel0Function = realApplyFunction;
			channel1Function = imagApplyFunction;
			fSequence.setChannelName(0, "Real");
			fSequence.setChannelName(1, "Imaginary");
		}
		
		AssignFunction2D assignFunction = null;
		if(!swap) //No Quadrant swapping
		{
			assignFunction = directAssign;
		}
		else //Swap quadrants
		{
			assignFunction = swapAssign;
		}

		for(int k = 0; k < _z; k++)
		{
			double[] fArray = new double[_w*_h*2];
			Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), 0, fArray, 0, _w*_h, sequence.isSignedDataType());
			
			// Computes 2D forward DFT of real data leaving the result in fArray
			// Because the result is stored in fArray, fArray must be of size rows*2*columns,
			// with only the first rows*columns elements filled with real data.
			fft.realForwardFull(fArray);

			IcyBufferedImage resultArray = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);
			double[][] resultData = resultArray.getDataXYCAsDouble();
					
			assignFunction.assign(fArray, resultData, _w, _h, channel0Function, channel1Function);

			resultArray.dataChanged();
			fSequence.setImage(0, k, resultArray);
		}

		return fSequence;
	}

	@Override
	public void clean() {
	}
}