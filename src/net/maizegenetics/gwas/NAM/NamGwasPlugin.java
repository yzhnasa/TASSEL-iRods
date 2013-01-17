package net.maizegenetics.gwas.NAM;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

public class NamGwasPlugin extends AbstractPlugin {
	private String mapFilename = null;
	private String residualFilename = null;
	private String rilmarkerFilename = null;
	private String founderFilename = null;
	private String modelFilename = null;
	private String stepFilename = null;
	private boolean randomizeSnpOrder = false;
	private boolean resample = true;
	private double enterLimit = 1e-6;
	private int iterations = 100;
	private int start = 1;
	private boolean threaded = false;
	
	private static final Logger myLogger = Logger.getLogger(NamGwasPlugin.class);

    public NamGwasPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

	@Override
	public DataSet performFunction(DataSet input) {
		if (mapFilename == null || residualFilename == null || rilmarkerFilename == null || founderFilename == null || modelFilename == null || stepFilename == null) {
			myLogger.info(getUsage());
			return null;
		}
		
		FileNames parameters = new FileNames();
		parameters.agpmap = new File(mapFilename);
		parameters.residuals = new File(residualFilename);
		parameters.namMarkersByChr = new File(rilmarkerFilename);
		parameters.snps = new File(founderFilename);
		parameters.chrmodel = new File(modelFilename);
		parameters.chrsteps = new File(stepFilename);
		parameters.enterlimit = enterLimit;
		parameters.iterations = iterations;
		parameters.startIteration = start;
		parameters.randomizeSnpOrder = randomizeSnpOrder;

		try {
			BufferedReader br = new BufferedReader(new FileReader(parameters.snps));
			br.readLine();
			String[] data = br.readLine().split("\t");
			parameters.chromosome = Integer.parseInt(data[2]);
			br.close();
		} catch(IOException e) {
			myLogger.error(e);
		}
		
		return null;
	}

    public void setMapFilename(String mapFilename) {
		this.mapFilename = mapFilename;
	}

	public void setResidualFilename(String residualFilename) {
		this.residualFilename = residualFilename;
	}

	public void setRilmarkerFilename(String rilmarkerFilename) {
		this.rilmarkerFilename = rilmarkerFilename;
	}

	public void setFounderFilename(String founderFilename) {
		this.founderFilename = founderFilename;
	}

	public void setRandomizeSnpOrder(boolean randomizeSnpOrder) {
		this.randomizeSnpOrder = randomizeSnpOrder;
	}

	public void setEnterLimit(double enterLimit) {
		this.enterLimit = enterLimit;
	}

	public void setResample(boolean resample) {
		this.resample = resample;
	}

	public void setIterations(int iterations) {
		this.iterations = iterations;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public void setThreaded(boolean threaded) {
		this.threaded = threaded;
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return null;
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equals("-c") || args[i].equalsIgnoreCase("-map")) {
				mapFilename = args[++i];
			}
			else if (args[i].equals("-t") || args[i].equalsIgnoreCase("-trait")) {
				residualFilename = args[++i];
			}
			else if (args[i].equals("-r") || args[i].equalsIgnoreCase("-rils")) {
				rilmarkerFilename = args[++i];
			}
			else if (args[i].equals("-f") || args[i].equalsIgnoreCase("-founders")) {
				founderFilename = args[++i];
			}
			else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-model")) {
				modelFilename = args[++i];
			}
			else if (args[i].equals("-s") || args[i].equalsIgnoreCase("-steps")) {
				stepFilename = args[++i];
			}
			else if (args[i].equals("-a") || args[i].equalsIgnoreCase("-randomize")) {
				if (args[++i].toUpperCase().startsWith("T")) randomizeSnpOrder = true;
			}
			else if (args[i].equals("-e") || args[i].equalsIgnoreCase("-enterlimit")) {
				enterLimit = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-i") || args[i].equalsIgnoreCase("-iterations")) {
				iterations = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-start")) {
				start = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-noresample")) {
				resample = false;
			}
			else if (args[i].equals("?")) myLogger.info(getUsage());
			else {
				myLogger.info(getUsage());
				return;
			}
		}
	}
	
	private String getUsage() {
		StringBuilder usage = new StringBuilder("The CallParentAllelesPlugin requires the following parameter:\n");
		usage.append("-c or -map : file name of RIL markers with map coordinates");
		usage.append("-t or -trait : file name of chromosome residuals for a trait");
		usage.append("-r or -rils : file name of the ril markers for this chromosome");
		usage.append("-f or -founders : file name of markers for this chromosome genotyped on founders to be projected on RILs");
		usage.append("-m or -model : file name of markers for this chromosome genotyped on founders to be projected on RILs");
		usage.append("-s or -steps : file name of markers for this chromosome genotyped on founders to be projected on RILs");
		usage.append("-a or -randomize : true if snps should be tested in random order (default = false)");
		usage.append("-e or -enterlimit : the largest p-value for which a new term will be added to the model (default = 1e-6");
		usage.append("-i or -iterations : the number of resampling iterations (default = 100");
		usage.append("-d or -start : the number of the first iteration in this sequence (default = 1");
		usage.append("-noresample : do not resample (does not take a parameter)");
		usage.append("-enablethreads : have application use multiple cores if available.");
		usage.append("? : print the parameter list.\n");

		return usage.toString();
	}
	
	 
}
