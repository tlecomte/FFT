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



	private Sequence FFT_3D(Sequence sequence, String swap, String display) {
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();
		int wc = (int) Math.ceil(_w/2);
		int hc = (int) Math.ceil(_h/2);
		int zc = (int) Math.ceil(_z/2);

		double[] fArray;
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

		if(swap == "No")
		{ //No Quadrant swapping. Leave as it is.
			for(int k = 0; k < _z; k++)
			{	
				IcyBufferedImage resultMatrix = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);			
				resultMatrix.setDataXY(0, Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), sequence.isSignedDataType()));//set buffered image to sequence 
				fSequence.setImage(0, k, resultMatrix);			
			}						
			fArray = fSequence.getDataCopyCXYZAsDouble(0);
			fft.complexForward(fArray);//Does only on half the data. To get the full transform use realForwardFull

			if(display=="Magnitude/Phase Pair")
			{
				for(int k = 0; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);


					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < _w; x++)
						{
							for(int y = 0; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[(x + (y * _w) + (k * _w * _h))*2 + 0], 2)+Math.pow(fArray[(x + (y * _w) + (k * _w * _h))*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[(x + (y * _w) + (k * _w * _h))*2 + 1], fArray[(x + (y * _w) + (k * _w * _h))*2 + 0]));
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
			}
			else
			{
				for(int k = 0; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);
					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < _w; x++)
						{
							for(int y = 0; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[(x + (y * _w) + (k * _w * _h))*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[(x + (y * _w) + (k * _w * _h))*2 + 1]);
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
			}
		}
		else
		{//Swap Quadrants
			for(int k = 0; k < _z; k++)
			{	
				IcyBufferedImage resultMatrix = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);			
				resultMatrix.setDataXY(0, Array1DUtil.arrayToDoubleArray(sequence.getDataXY(0, k, 0), sequence.isSignedDataType()));//set buffered image to sequence 
				fSequence.setImage(0, k, resultMatrix);			
			}						
			fArray = fSequence.getDataCopyCXYZAsDouble(0);
			fft.complexForward(fArray);//Does only on half the data. To get the full transform use realForwardFull

			if(display=="Magnitude/Phase Pair")
			{
				for(int k = 0; k < zc+1; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);


					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1], fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0]));
							}
						}
					}

					finally
					{
						resultMatrix.endUpdate();
					}
				}
				for(int k = zc+1; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);


					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], fArray[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], fArray[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]));
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], fArray[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]));
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0], 2)+Math.pow(fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], 2)));
								resultMatrix.setDataAsDouble(x, y, 1, Math.atan2(fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1], fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]));
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
			}
			else
			{
				for(int k = 0; k < zc+1; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);
					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((_w+(wc-x)) + (hc-y) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (zc-k) * _w * _h)*2 + 1]);
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}
				for(int k = zc+1; k < _z; k++)
				{			
					IcyBufferedImage resultMatrix = fSequence.getImage(0, k);
					resultMatrix.beginUpdate();
					try
					{
						for(int x = 0; x < wc+1; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((wc-x) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1]);
							}
						}
						for(int x = wc+1; x < _w; x++)
						{
							for(int y = 0; y < hc+1; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((_w+(wc-x)) + (hc-y) * _w + (_z+(zc-k)) * _w * _h)*2 + 1]);
							}
							for(int y = hc+1; y < _h; y++)
							{
								resultMatrix.setDataAsDouble(x, y, 0, fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 0]);
								resultMatrix.setDataAsDouble(x, y, 1, fArray[((_w+(wc-x)) + (_h+(hc-y)) * _w + (_z+(zc-k)) * _w * _h)*2 + 1]);
							}
						}
					}
					finally
					{
						resultMatrix.endUpdate();
					}
				}

			}

		}
		addSequence(fSequence);
		return fSequence;
	}
	
	interface ApplyFunction {
		double apply0(double real, double imag);
		double apply1(double real, double imag);
	}
	
	ApplyFunction realImageApplyFunction = new ApplyFunction()
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
	
	interface AssignFunction {
		void assign(double[] in, double[][] out, int _w, int _h, ApplyFunction function);
	}
	
	AssignFunction directAssign = new AssignFunction()
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
	
	AssignFunction swapAssign = new AssignFunction() {
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
			
			ApplyFunction applyFunction = null;
			if(display=="Magnitude/Phase Pair")
			{
				applyFunction = magnitudeAngleApplyFunction;
			}
			else // Real/Imaginary Pair
			{
				applyFunction = realImageApplyFunction;
			}
			
			
			AssignFunction assignFunction = null;
			if(swap == "No") //No Quadrant swapping
			{
				assignFunction = directAssign;
			}
			else //Swap quadrants
			{
				assignFunction = swapAssign;
			}
			
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