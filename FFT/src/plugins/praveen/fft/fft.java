package plugins.praveen.fft;

import cern.colt.function.tdouble.DoubleDoubleFunction;
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

public class FFT extends EzPlug implements Block {

	enum FFTDims {
		FFT_2D("2D (xy)"), FFT_3D("3D (xyz)");
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
		{		
			fSequence = FFT_2D(sequence, swap.getValue(), outputType.getValue());
		}
		else
		{
			if (sequence.getSizeZ() >= 2)
			{
				fSequence = FFT_3D(sequence, swap.getValue(), outputType.getValue());
			}
			else
			{
				System.err.println("Sequence depth is 1, so computing 2D FFT instead of 3D.");
				fSequence = FFT_2D(sequence, swap.getValue(), outputType.getValue());
			}
		}
		
		if (!isHeadLess()) {
			addSequence(fSequence);
		}
		
		fSequenceVar.setValue(fSequence);
	}

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

		DoubleDoubleFunction channel0ApplyFunction = null;
		DoubleDoubleFunction channel1ApplyFunction = null;
		if(outputType == FFTOutputType.MAGNITUDE_PHASE)
		{
			channel0ApplyFunction = new ComplexFunctions.Magnitude();
			channel1ApplyFunction = new ComplexFunctions.Angle();
			fSequence.setChannelName(0, "Magnitude");
			fSequence.setChannelName(1, "Phase");
		}
		else
		{
			channel0ApplyFunction = new ComplexFunctions.Real();
			channel1ApplyFunction = new ComplexFunctions.Imag();
			fSequence.setChannelName(0, "Real");
			fSequence.setChannelName(1, "Imaginary");
		}
		
		AssignFunction3D assignFunction = null;
		if(!swap)  // No Quadrant swapping. Leave as it is.
		{
			assignFunction = new AssignFunctions.DirectAssign3D();
		}
		else
		{
			assignFunction = new AssignFunctions.SwapAssign3D(); // Swap Quadrants
		}
		
		assignFunction.assign(fArray, resultData, _w, _h, _z, 0, channel0ApplyFunction);
		assignFunction.assign(fArray, resultData, _w, _h, _z, 1, channel1ApplyFunction);

		fSequence.dataChanged();
		
		return fSequence;
	}

	private Sequence FFT_2D(Sequence sequence, boolean swap, FFTOutputType outputType) 
	{
		Sequence fSequence = new Sequence();
		fSequence.setName("Fourier Transform 2D");
		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();

		final DoubleFFT_2D fft = new DoubleFFT_2D(_h, _w);
		
		DoubleDoubleFunction channel0Function = null;
		DoubleDoubleFunction channel1Function = null;
		if(outputType == FFTOutputType.MAGNITUDE_PHASE)
		{
			channel0Function = new ComplexFunctions.Magnitude();
			channel1Function = new ComplexFunctions.Angle();
			fSequence.setChannelName(0, "Magnitude");
			fSequence.setChannelName(1, "Phase");
		}
		else // Real/Imaginary Pair
		{
			channel0Function = new ComplexFunctions.Real();
			channel1Function = new ComplexFunctions.Imag();
			fSequence.setChannelName(0, "Real");
			fSequence.setChannelName(1, "Imaginary");
		}
		
		AssignFunction2D assignFunction = null;
		if(!swap) //No Quadrant swapping
		{
			assignFunction = new AssignFunctions.DirectAssign2D();
		}
		else //Swap quadrants
		{
			assignFunction = new AssignFunctions.SwapAssign2D();
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
					
			assignFunction.assign(fArray, resultData[0], _w, _h, channel0Function);
			assignFunction.assign(fArray, resultData[1], _w, _h, channel1Function);

			resultArray.dataChanged();
			fSequence.setImage(0, k, resultArray);
		}

		return fSequence;
	}

	@Override
	public void clean() {
	}
}