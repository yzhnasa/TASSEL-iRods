package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.pal.alignment.TBitAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

public class ViterbiAlgorithmPlugin extends AbstractPlugin {

	private static final Logger myLogger = Logger.getLogger(CallParentAllelesPlugin.class);
	private boolean fillGapsInAlignment = true;
	private double probHeterozygous = 0.07;
	
	public ViterbiAlgorithmPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> theData = input.getDataOfType(PopulationData.class);
		for (Datum data:theData) {
			PopulationData family = (PopulationData) data.getData();
			TBitAlignment tba;
			if (family.imputed instanceof TBitAlignment) {
				tba = (TBitAlignment) family.imputed;
			} else {
				tba = TBitAlignment.getInstance(family.imputed);
			}
			
			double phet;
			if (family.inbredCoef >= 0 && family.inbredCoef <= 1) {
				phet = (1 - family.inbredCoef) / 2;
			} else {
				phet = probHeterozygous;
			}
			
			family.imputed = NucleotideImputationUtils.imputeUsingViterbiFiveState(tba, phet, family.name);
			
			if (fillGapsInAlignment) NucleotideImputationUtils.fillGapsInAlignment(family);
		}
		
		DataSet resultDS = new DataSet(theData, this);
		fireDataSetReturned(new PluginEvent(resultDS, ViterbiAlgorithmPlugin.class));
		return resultDS;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
//			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equals("-g") || args[i].equalsIgnoreCase("-fillgaps")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) fillGapsInAlignment = true;
				else fillGapsInAlignment = false;
			}
			else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-phet")) {
				probHeterozygous = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("?")) myLogger.error(getUsage());
		}
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "Viterbi";
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	public void setFillGapsInAlignment(boolean fillGapsInAlignment) {
		this.fillGapsInAlignment = fillGapsInAlignment;
	}

	public void setProbHeterozygous(double probHeterozygous) {
		this.probHeterozygous = probHeterozygous;
	}

	private String getUsage() {
		StringBuilder usage = new StringBuilder("The ViterbiAlgorithmPlugin can take the following optional parameters:\n");
		usage.append("-g or -fillgaps : true if missing values between SNPs from the same parent should be imputed to that parent, false otherwise (default = true)\n");
		usage.append("-h or -phet : expected frequency of heterozygous loci (default = 0.07). If the inbreeding coefficient is specified in the pedigree file that will used to calculate this value.\n");
		usage.append("? : print the parameter list.\n");

		return usage.toString();
	}

}
