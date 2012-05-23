package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.io.File;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

public class WritePopulationAlignmentPlugin extends AbstractPlugin {
	private static final Logger myLogger = Logger.getLogger(WritePopulationAlignmentPlugin.class);
	boolean mergeAlignments = true;
	boolean writeParentCalls = true;
	boolean outputDiploid = false;
	String baseFileName;
	
	public WritePopulationAlignmentPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> theData = input.getDataOfType(PopulationData.class);
		if (mergeAlignments) {
			String filename = baseFileName + ".hmp.txt";
			Alignment[] allOfTheAlignments = new Alignment[theData.size()];
			int count = 0;
			for (Datum data:theData) {
				if (writeParentCalls) allOfTheAlignments[count++] = ((PopulationData) data.getData()).imputed;
				else allOfTheAlignments[count++] = ((PopulationData) data.getData()).original;
			}
			Alignment alignment = MutableSingleEncodeAlignment.getInstance(allOfTheAlignments);
			ExportUtils.writeToHapmap(alignment, outputDiploid, filename, '\t', null);
		} else {
			for (Datum datum:theData) {
				PopulationData family = (PopulationData) datum.getData();
				String filename = baseFileName + ".family." + family.name + ".hmp.txt";
				if (writeParentCalls) ExportUtils.writeToHapmap(family.imputed, outputDiploid, filename, '\t', null);
				else ExportUtils.writeToHapmap(family.original, outputDiploid, filename, '\t', null);
			}
		}
		
		fireDataSetReturned(new PluginEvent(input, input.getCreator().getClass()));
		return input;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equals("-f") || args[i].equalsIgnoreCase("-file")) {
				baseFileName = args[++i];
			}
			else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-merge")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) mergeAlignments = true;
				else mergeAlignments = false;
			}
			else if (args[i].equals("-p") || args[i].equalsIgnoreCase("-parentCalls")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) writeParentCalls = true;
				else writeParentCalls = false;
			}
			else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-diploid")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) writeParentCalls = true;
				else writeParentCalls = false;
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
		return "Write Populations";
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	private String getUsage() {
		StringBuilder usage = new StringBuilder("The WritePopulationAlignmentPlugin requires the following parameter:\n");
		usage.append("-f or -file : The base file name for the ouput. .hmp.txt will be appended.\n");
		usage.append("The following parameters are optional:\n");
		usage.append("-m or -merge : if true families are merged into a single file, if false each family is output to a separate file (default = true)");
		usage.append("-p or -parentCalls : if true output is A/C/M where parent 1 is A and parent 2 is M, if false outputs nucleotides (default = true)");
		usage.append("-d or -diploid : if true output is AA/CC/AC, if false output is A/C/M");
		usage.append("? : print the parameter list.");

		return usage.toString();
	}

}
