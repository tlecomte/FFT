package plugins.praveen.fft;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;

public class fft extends EzPlug {

	EzVarSequence input = new EzVarSequence("Input");
	EzVarText	ndims = new EzVarText("Type", new String[] { "2D", "3D" }, 0, false);
	EzVarText	display = new EzVarText("Display as", new String[] {  "Magnitude/Phase Pair", "Real/Imaginary Pair" }, 0, false);
	EzVarText	swap = new EzVarText("Swap Quadrants?", new String[] { "Yes", "No" }, 1, false);

	@Override
	protected void initialize() {
		super.addEzComponent(input);
		super.addEzComponent(ndims);
		super.addEzComponent(swap);
		super.addEzComponent(display);		
		super.setTimeDisplay(true);
	}

	@Override
	protected void execute() {
		Sequence sequence = input.getValue();

		if(ndims.getValue()=="2D")		
			FFT_2D(sequence, swap.getValue(), display.getValue());	
		else
			FFT_3D(sequence, swap.getValue(), display.getValue());
		//MessageDialog.showDialog("FFT3D not implemented yet !");	
	}
	
	interface ApplyFunction {
		double apply0(double real, double imag);
		double apply1(double real, double imag);
	}
	
	ApplyFunction realImagApplyFunction = new ApplyFunction()
	{
		public double apply0(double real, double imag)
		{
			return real;
		}
		public double apply1(double real, double imag)
		{
			return imag;
		}
	};
	
	ApplyFunction magnitudeAngleApplyFunction = new ApplyFunction()
	{
		public double apply0(double real, double imag)
		{
			return Math.sqrt(Math.pow(real, 2) + Math.pow(imag, 2));
		}
		public double apply1(double real, double imag)
		{
			return Math.atan2(imag, real);
		}
	};

	interface AssignFunction3D {
		void assign(double[] in, double[][][] out, int _w, int _h, int _z, ApplyFunction function);
	}
	
	// function that walks the 3D FFT from JTransforms and fills the sequence data array
	// this is the version that does not swap the quadrants
	AssignFunction3D directAssign3D = new AssignFunction3D()
	{
		public void assign(double[] in, double[][][] out, int _w, int _h, int _z, ApplyFunction function) {
			for(int k = 0; k < _z; k++)
			{			
				for(int x = 0; x < _w; x++)
				{
					for(int y = 0; y < _h; y++)
					{
						double real = in[(x + (y * _w) + (k * _w * _h))*2 + 0];
						double imag = in[(x + (y * _w) + (k * _w * _h))*2 + 1];					
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
				}
			}
		}	
	};

	// function that walks the 3D FFT from JTransforms and fills the sequence data array
	// this is the version that swaps the quadrants
	AssignFunction3D swapAssign3D = new AssignFunction3D()
	{
		public void assign(double[] in, double[][][] out, int _w, int _h, int _z, ApplyFunction function) {
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
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0];
						double imag = in[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
				}
				for(int x = wc+1; x < _w; x++)
				{
					for(int y = 0; y < hc+1; y++)
					{
						double real = in[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
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
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
						double imag = in[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
				}
				for(int x = wc+1; x < _w; x++)
				{
					for(int y = 0; y < hc+1; y++)
					{
						double real = in[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0];
						double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1];
						out[k][0][x + _w*y] = function.apply0(real, imag);
						out[k][1][x + _w*y] = function.apply1(real, imag);
					}
				}
			}
		}
	};

	private Sequence FFT_3D(Sequence sequence, String swap, String display) {
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();

		final DoubleFFT_3D fft = new DoubleFFT_3D(_z, _h, _w);
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 3D");

		if(display=="Magnitude/Phase Pair")
		{
			fSequence.setChannelName(0, "Magnitude");
			fSequence.setChannelName(1, "Phase");
		}
		else // Real/Imaginary Pair
		{
			fSequence.setChannelName(0, "Real");
			fSequence.setChannelName(1, "Imaginary");
		}

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

		ApplyFunction applyFunction = null;
		if(display=="Magnitude/Phase Pair")
		{
			applyFunction = magnitudeAngleApplyFunction;
		}
		else
		{
			applyFunction = realImagApplyFunction;
		}
		
		AssignFunction3D assignFunction = null;
		if(swap == "No")  // No Quadrant swapping. Leave as it is.
		{
			assignFunction = directAssign3D;
		}
		else
		{
			assignFunction = swapAssign3D; // Swap Quadrants
		}
		
		assignFunction.assign(fArray, resultData, _w, _h, _z, applyFunction);

		fSequence.dataChanged();
		
		addSequence(fSequence);
		return fSequence;
	}
	
	interface Assign2DFunction {
		void assign(double[] in, double[][] out, int _w, int _h, ApplyFunction function);
	}
	
	// function that walks the 2D FFT from JTransforms and fills the sequence data array
	// this is the version that does not swap the quadrants
	Assign2DFunction directAssign = new Assign2DFunction()
	{
		public void assign(double[] in, double[][] out, int _w, int _h, ApplyFunction function) {
			for (int i = 0; i < in.length/2; i++)
			{
				double real = in[2*i];
				double imag = in[2*i + 1];
				
				out[0][i] = function.apply0(real, imag);
				out[1][i] = function.apply1(real, imag);
			}	
		}	
	};
	
	// function that walks the 2D FFT from JTransforms and fills the sequence data array
	// this is the version that swaps the quadrants
	Assign2DFunction swapAssign = new Assign2DFunction() {
		public void assign(double[] in, double[][] out, int _w, int _h, ApplyFunction function)
		{
			int wc = (int) Math.ceil(_w/2);
			int hc = (int) Math.ceil(_h/2);
			
			for(int x = 0; x < (wc+1); x++)
			{
				for(int y = 0; y < (hc+1); y++)
				{
					double real = in[((wc-x) + (hc-y) * _w)*2 + 0];
					double imag = in[((wc-x) + (hc-y) * _w)*2 + 1];
					
					out[0][x + _w*y] = function.apply0(real, imag);
					out[1][x + _w*y] = function.apply1(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((wc-x) + (_h+(hc-y)) * _w)*2 + 0];
					double imag = in[((wc-x) + (_h+(hc-y)) * _w)*2 + 1];
					
					out[0][x + _w*y] = function.apply0(real, imag);
					out[1][x + _w*y] = function.apply1(real, imag);
				}
		
			}
			for(int x = (wc+1); x < _w; x++)
			{
				for(int y = 0; y < (hc+1); y++)
				{
					double real = in[((_w+(wc-x)) + (hc-y) * _w)*2 + 0];
					double imag = in[((_w+(wc-x)) + (hc-y) * _w)*2 + 1];
					
					out[0][x + _w*y] = function.apply0(real, imag);
					out[1][x + _w*y] = function.apply1(real, imag);
				}
				for(int y = hc+1; y < _h; y++)
				{
					double real = in[((_w+(wc-x)) + (_h+(hc-y)) * _w)*2 + 0];
					double imag = in[((_w+(wc-x)) + (_h+(hc-y)) * _w)*2 + 1];
					
					out[0][x + _w*y] = function.apply0(real, imag);
					out[1][x + _w*y] = function.apply1(real, imag);
				}
			}
		}
	};

	private Sequence FFT_2D(Sequence sequence, String swap, String display) 
	{
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 2D");
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();
		
		if(display=="Magnitude/Phase Pair")
		{
			fSequence.setChannelName(0, "Magnitude");
			fSequence.setChannelName(1, "Phase");
		}
		else // Real/Imaginary Pair
		{
			fSequence.setChannelName(0, "Real");
			fSequence.setChannelName(1, "Imaginary");
		}

		final DoubleFFT_2D fft = new DoubleFFT_2D(_h, _w);
		
		ApplyFunction applyFunction = null;
		if(display=="Magnitude/Phase Pair")
		{
			applyFunction = magnitudeAngleApplyFunction;
		}
		else // Real/Imaginary Pair
		{
			applyFunction = realImagApplyFunction;
		}
		
		
		Assign2DFunction assignFunction = null;
		if(swap == "No") //No Quadrant swapping
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
					
			assignFunction.assign(fArray, resultData, _w, _h, applyFunction);

			resultArray.dataChanged();
			fSequence.setImage(0, k, resultArray);
		}

		addSequence(fSequence);

		return fSequence;
	}

	@Override
	public void clean() {
	}
}