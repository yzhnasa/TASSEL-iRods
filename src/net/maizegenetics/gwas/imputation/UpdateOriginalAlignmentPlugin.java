package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.pal.alignment.TBitAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class UpdateOriginalAlignmentPlugin extends AbstractPlugin {
	//this plugin does not take any user-defined parameters
	
	public UpdateOriginalAlignmentPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> theData = input.getDataOfType(PopulationData.class);
		for (Datum data:theData) {
			PopulationData family = (PopulationData) data.getData();
			NucleotideImputationUtils.updateSnpAlignment(family);
		}
		return new DataSet(theData, this);
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "ImputeOriginal";
	}

	@Override
	public String getToolTipText() {
		return null;
	}

}
