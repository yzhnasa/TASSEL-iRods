/**************************** Zhong Yang *********************************/
// This class is for implement "Load iRods" item in Data menu.
package net.maizegenetics.analysis.data;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JRootPane;
import javax.swing.tree.TreePath;

import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.dna.map.TOPMUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.ReadPolymorphismUtils;
import net.maizegenetics.dna.snp.ReadSequenceAlignmentUtils;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.tassel.TASSELMainFrame;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.trait.ReadPhenotypeUtils;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Report;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;
import org.bio5.irods.iplugin.bean.TapasCoreFunctions;
import org.bio5.irods.iplugin.utilities.IrodsUtilities;
import org.bio5.irods.iplugin.views.IPlugin_OpenImage;
import org.bio5.irods.iplugin.views.MainWindow;
import org.irods.jargon.core.exception.JargonException;
import org.irods.jargon.core.pub.io.IRODSFile;
import org.irods.jargon.core.pub.io.IRODSFileFactory;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

public class IRodsFileLoadPlugin extends AbstractPlugin {

	private PlinkLoadPlugin myPlinkLoadPlugin = null;
	private ProjectionLoadPlugin myProjectionLoadPlugin = null;
	private String[] myOpenFiles;
	private ArrayList myOpenFilesList = new ArrayList<String>();
	private TasselFileType myFileType = TasselFileType.Unknown;
	private static final Logger myLogger = Logger
			.getLogger(FileLoadPlugin.class);
	private MainWindow irodsFileChooser = null;
	private IRODSFileFactory iRODSFileFactory;
	IRODSFile sourceIrodsFilePath = null;
	private List result = new ArrayList();
	private ArrayList<String> alreadyLoaded = new ArrayList();

	public IRodsFileLoadPlugin(Frame parentFame, boolean isInterative) {
		super(parentFame, isInterative);
	}
	
	public IRodsFileLoadPlugin(Frame parentFame, boolean isInterative,
			PlinkLoadPlugin plinkLoadPlugin,
			ProjectionLoadPlugin projectionLoadPlugin) {
		super(parentFame, isInterative);
		myPlinkLoadPlugin = plinkLoadPlugin;
		myProjectionLoadPlugin = projectionLoadPlugin;
	}

	public IRodsFileLoadPlugin(Frame parentFame, boolean isInterative,
			PlinkLoadPlugin plinkLoadPlugin,
			ProjectionLoadPlugin projectionLoadPlugin,
			MainWindow irodsFileChooser) {
		super(parentFame, isInterative);
		myPlinkLoadPlugin = plinkLoadPlugin;
		myProjectionLoadPlugin = projectionLoadPlugin;
		this.irodsFileChooser = irodsFileChooser;
	}
	
	public void setIRodsFileChooser(MainWindow irodsFileChooser){
		this.irodsFileChooser = irodsFileChooser;
	}
	
	public TasselFileType getFileType(){
		return myFileType;
	}
	
	public Logger getLogger(){
		return myLogger;
	}
	
	public PlinkLoadPlugin getPlinkLoadPlugin(){
		return myPlinkLoadPlugin;
	}
	
	public void setResult(DataSet tds){
		result.add(tds);
	}
	
	public void setAlreadyLoaded(String myOpenFile){
		alreadyLoaded.add(myOpenFile);
	}

	public DataSet performFunction(DataSet input) {

		try {

			if (isInteractive()) {
				FileLoadPluginDialog theDialog = new FileLoadPluginDialog();
				theDialog.setLocationRelativeTo(getParentFrame());
				theDialog.setVisible(true);
				if (theDialog.isCancel()) {
					return null;
				}
				myFileType = theDialog.getTasselFileType();

				if (myFileType == TasselFileType.Plink) {
					return myPlinkLoadPlugin.performFunction(null);
				}

				if (myFileType == TasselFileType.ProjectionAlignment) {
					return myProjectionLoadPlugin.performFunction(null);
				}
				getOpenFilesByChooser();
				theDialog.dispose();
			}

			return DataSet.getDataSet(result, this);

		} finally {
			fireProgress(100);
		}

	}

	public DataSet guessAtUnknowns(String filename) {

		TasselFileType guess = TasselFileType.Sequence;
		DataSet tds = null;

		try {
			BufferedReader br = null;
			if (filename.startsWith("http")) {
				URL url = new URL(filename);
				br = new BufferedReader(new InputStreamReader(url.openStream()));
			} else {
				br = new BufferedReader(new FileReader(filename));
			}

			String line1 = br.readLine().trim();
			String[] sval1 = line1.split("\\s");
			String line2 = br.readLine().trim();
			String[] sval2 = line2.split("\\s");
			if (line1.startsWith("<") || line1.startsWith("#")) {
				boolean isTrait = false;
				boolean isMarker = false;
				boolean isNumeric = false;
				boolean isMap = false;
				Pattern tagPattern = Pattern.compile("[<>\\s]+");
				String[] info1 = tagPattern.split(line1);
				String[] info2 = tagPattern.split(line2);
				if (info1.length > 1) {
					if (info1[1].toUpperCase().startsWith("MARKER")) {
						isMarker = true;
					} else if (info1[1].toUpperCase().startsWith("TRAIT")) {
						isTrait = true;
					} else if (info1[1].toUpperCase().startsWith("NUMER")) {
						isNumeric = true;
					} else if (info1[1].toUpperCase().startsWith("MAP")) {
						isMap = true;
					}
				}
				if (info2.length > 1) {
					if (info2[1].toUpperCase().startsWith("MARKER")) {
						isMarker = true;
					} else if (info2[1].toUpperCase().startsWith("TRAIT")) {
						isTrait = true;
					} else if (info2[1].toUpperCase().startsWith("NUMER")) {
						isNumeric = true;
					} else if (info2[1].toUpperCase().startsWith("MAP")) {
						isMap = true;
					}
				} else {
					guess = null;
					String inline = br.readLine();
					while (guess == null
							&& inline != null
							&& (inline.startsWith("#") || inline
									.startsWith("<"))) {
						if (inline.startsWith("<")) {
							String[] info = tagPattern.split(inline);
							if (info[1].toUpperCase().startsWith("MARKER")) {
								isMarker = true;
							} else if (info[1].toUpperCase()
									.startsWith("TRAIT")) {
								isTrait = true;
							} else if (info[1].toUpperCase()
									.startsWith("NUMER")) {
								isNumeric = true;
							} else if (info[1].toUpperCase().startsWith("MAP")) {
								isMap = true;
							}
						}
					}
				}
				if (isTrait || (isMarker && isNumeric)) {
					guess = TasselFileType.Phenotype;
				} else if (isMap) {
					guess = TasselFileType.GeneticMap;
				} else {
					throw new IOException(
							"Improperly formatted header. Data will not be imported.");
				}
			} else if ((line1.startsWith(">")) || (line1.startsWith(";"))) {
				guess = TasselFileType.Fasta;
			} else if (sval1.length == 1) {
				guess = TasselFileType.SqrMatrix;
			} else if ((line1.startsWith("#Nexus"))
					|| (line1.startsWith("#NEXUS"))
					|| (line1.startsWith("CLUSTAL"))
					|| ((sval1.length == 2) && (sval2.length == 2))) {
				guess = TasselFileType.Sequence;
			}

			myLogger.info("guessAtUnknowns: type: " + guess);
			tds = processDatum(filename, guess);

			br.close();
		} catch (Exception e) {
		}

		return tds;

	}

	public DataSet processDatum(String inFile, TasselFileType theFT) {
		Object result = null;
		String suffix = null;
		try {
			switch (theFT) {
			case Hapmap: {
				suffix = FileLoadPlugin.FILE_EXT_HAPMAP;
				if (inFile.endsWith(".gz")) {
					suffix = FileLoadPlugin.FILE_EXT_HAPMAP_GZ;
				}
				result = ImportUtils.readFromHapmap(inFile, this);
				break;
			}
			case HDF5: {
				IHDF5Reader reader = HDF5Factory.openForReading(inFile);
				boolean t4HDF5 = HDF5Utils.isTASSEL4HDF5Format(HDF5Factory
						.openForReading(inFile));
				reader.close();
				if (t4HDF5) {
					DialogUtils
							.showError(
									"This file is TASSEL4 HDF5 format. It will be converted to TASSEL5 "
											+ "HDF5 format with the t5.h5 suffix added.  This may take a few minutes.",
									getParentFrame());
					String newInfile = inFile.replace(".h5", ".t5.h5");
					MigrateHDF5FromT4T5.copyGenotypes(inFile, newInfile);
					inFile = newInfile;
				}
				suffix = FileLoadPlugin.FILE_EXT_HDF5;
				result = ImportUtils.readGuessFormat(inFile);
				break;
			}
			case VCF: {
				suffix = FileLoadPlugin.FILE_EXT_VCF;
				if (inFile.endsWith(".gz")) {
					suffix = FileLoadPlugin.FILE_EXT_VCF + ".gz";
				}
				result = ImportUtils.readFromVCF(inFile, this);
				break;
			}
			case Sequence: {
				result = ReadSequenceAlignmentUtils.readBasicAlignments(inFile,
						40);
				break;
			}
			case Fasta: {
				result = ImportUtils.readFasta(inFile);
				break;
			}
			case SqrMatrix: {
				result = ReadDistanceMatrix.readDistanceMatrix(inFile);
				break;
			}
			case Phenotype: {
				result = ReadPhenotypeUtils.readGenericFile(inFile);
				break;
			}
			case GeneticMap: {
				result = ReadPolymorphismUtils.readGeneticMapFile(inFile);
				break;
			}
			case Table: {
				result = TableReportUtils
						.readDelimitedTableReport(inFile, "\t");
				break;
			}
			case TOPM: {
				result = TOPMUtils.readTOPM(inFile);
				break;
			}
			}
		} catch (Exception e) {

			e.printStackTrace();
			StringBuilder builder = new StringBuilder();
			builder.append("Error loading: ");
			builder.append(inFile);
			builder.append("\n");
			builder.append(Utils.shortenStrLineLen(
					ExceptionUtils.getExceptionCauses(e), 50));
			String str = builder.toString();
			if (isInteractive()) {
				DialogUtils.showError(str, getParentFrame());
			} else {
				myLogger.error(str);
			}

		}
		if (result != null) {
			String theComment = "";
			if (result instanceof Report) {
				StringWriter sw = new StringWriter();
				((Report) result).report(new PrintWriter(sw));
				theComment = sw.toString();
			}

			String name = Utils.getFilename(inFile, suffix);

			Datum td = new Datum(name, result, theComment);
			// todo need to add logic of directories.
			DataSet tds = new DataSet(td, this);
			return tds;
		}
		return null;
	}
	
	private void getOpenFilesByChooser() {
		if(IPlugin_OpenImage.getIrodsImagej().getMainWindowStatus() == true){
			irodsFileChooser = IPlugin_OpenImage.getIrodsImagej().getMainWindow();
			irodsFileChooser.setVisible(true);
		}else{
			IPlugin_OpenImage.getIrodsImagej().getMainWindow().setVisible(true);
			IPlugin_OpenImage.getIrodsImagej().getMainWindow().setVisibleFlag(true);
			if(IPlugin_OpenImage.getIrodsImagej().getMainWindowStatus() == true)
				IPlugin_OpenImage.getIrodsImagej().getTASSELMainFrame().setLoadDataFromIRodsMenuItemEnable();
		}
	}
	
	public String[] getOpenFiles() {
		return myOpenFiles;
	}

	public void setOpenFiles(File[] openFiles) {

		if ((openFiles == null) || (openFiles.length == 0)) {
			myOpenFiles = null;
			return;
		}

		myOpenFiles = new String[openFiles.length];
		for (int i = 0; i < openFiles.length; i++) {
			myOpenFiles[i] = openFiles[i].getPath();
		}

	}
	
	public void setOpenFiles(File openFile) {

		if (openFile == null) {
			myOpenFiles = null;
			return;
		}

		myOpenFilesList.add(openFile.getPath());
		myOpenFiles = (String[]) myOpenFilesList.toArray(new String[myOpenFilesList.size()]);
	}

	@Override
	public ImageIcon getIcon() {
		 URL imageURL = IRodsFileLoadPlugin.class.getResource("/net/maizegenetics/analysis/images/iplant.gif");
	        if (imageURL == null) {
	            return null;
	        } else {
	            return new ImageIcon(imageURL);
	        }
	}

	@Override
	public String getButtonName() {
		return "Load iRods";
	}

	@Override
	public String getToolTipText() {
		return "Load data from files on iRods server.";
	}

}

class FileLoadPluginDialog extends JDialog {

	boolean isCancel = true;
	ButtonGroup conversionButtonGroup = new ButtonGroup();
	JRadioButton hapMapRadioButton = new JRadioButton("Load Hapmap");
	JRadioButton hdf5RadioButton = new JRadioButton("Load HDF5");
	JRadioButton vcfRadioButton = new JRadioButton("Load VCF");
	JRadioButton plinkRadioButton = new JRadioButton("Load Plink");
	JRadioButton sequenceAlignRadioButton = new JRadioButton("Load Phylip");
	JRadioButton fastaRadioButton = new JRadioButton("Load FASTA File");
	JRadioButton numericalRadioButton = new JRadioButton(
			"Load Numerical (trait, covariates, or factors)");
	JRadioButton loadMatrixRadioButton = new JRadioButton(
			"Load Square Numerical Matrix (i.e. kinship)");
	JRadioButton guessRadioButton = new JRadioButton("Make Best Guess");
	JRadioButton projectionAlignmentRadioButton = new JRadioButton(
			"Load Projection Alignment");
	JRadioButton geneticMapRadioButton = new JRadioButton("Load a Genetic Map");
	JRadioButton tableReportRadioButton = new JRadioButton(
			"Load a Table Report");
	JRadioButton topmRadioButton = new JRadioButton(
			"Load a TOPM (Tags on Physical Map)");

	public FileLoadPluginDialog() {
		super((Frame) null, "File Loader", true);
		try {
			jbInit();
			pack();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception {

		setTitle("File Loader");
		setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
		setUndecorated(false);
		getRootPane().setWindowDecorationStyle(JRootPane.NONE);

		Container contentPane = getContentPane();

		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
		contentPane.setLayout(layout);

		JPanel main = getMain();

		contentPane.add(main);

		pack();

		setResizable(false);

		conversionButtonGroup.add(projectionAlignmentRadioButton);
		conversionButtonGroup.add(hapMapRadioButton);
		conversionButtonGroup.add(hdf5RadioButton);
		conversionButtonGroup.add(vcfRadioButton);
		conversionButtonGroup.add(plinkRadioButton);
		conversionButtonGroup.add(sequenceAlignRadioButton);
		conversionButtonGroup.add(fastaRadioButton);
		conversionButtonGroup.add(loadMatrixRadioButton);
		conversionButtonGroup.add(numericalRadioButton);
		conversionButtonGroup.add(geneticMapRadioButton);
		conversionButtonGroup.add(tableReportRadioButton);
		conversionButtonGroup.add(topmRadioButton);
		conversionButtonGroup.add(guessRadioButton);
		guessRadioButton.setSelected(true);

	}

	private JPanel getMain() {

		JPanel inputs = new JPanel();
		BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
		inputs.setLayout(layout);
		inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);

		inputs.add(Box.createRigidArea(new Dimension(1, 10)));

		inputs.add(getLabel());

		inputs.add(Box.createRigidArea(new Dimension(1, 10)));

		inputs.add(getOptionPanel());

		inputs.add(Box.createRigidArea(new Dimension(1, 10)));

		inputs.add(getButtons());

		inputs.add(Box.createRigidArea(new Dimension(1, 10)));

		return inputs;

	}

	private JPanel getLabel() {

		JPanel result = new JPanel();
		BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
		result.setLayout(layout);
		result.setAlignmentX(JPanel.CENTER_ALIGNMENT);

		JLabel jLabel1 = new JLabel("Choose File Type to Load.");
		jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
		result.add(jLabel1);

		return result;

	}

	private JPanel getOptionPanel() {

		JPanel result = new JPanel();
		BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
		result.setLayout(layout);
		result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
		result.setBorder(BorderFactory.createEtchedBorder());

		result.add(hapMapRadioButton);
		result.add(hdf5RadioButton);
		result.add(vcfRadioButton);
		result.add(plinkRadioButton);
		result.add(projectionAlignmentRadioButton);
		result.add(sequenceAlignRadioButton);
		result.add(fastaRadioButton);
		result.add(numericalRadioButton);
		result.add(loadMatrixRadioButton);
		result.add(geneticMapRadioButton);
		result.add(tableReportRadioButton);
		result.add(topmRadioButton);
		result.add(guessRadioButton);

		result.add(Box.createRigidArea(new Dimension(1, 20)));

		return result;

	}

	private JPanel getButtons() {

		JButton okButton = new JButton();
		JButton cancelButton = new JButton();

		cancelButton.setText("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				cancelButton_actionPerformed(e);
			}
		});

		okButton.setText("OK");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				okButton_actionPerformed(e);
			}
		});

		JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

		result.add(okButton);

		result.add(cancelButton);

		return result;

	}

	public FileLoadPlugin.TasselFileType getTasselFileType() {
		if (hapMapRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.Hapmap;
		}
		if (hdf5RadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.HDF5;
		}
		if (vcfRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.VCF;
		}
		if (plinkRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.Plink;
		}
		if (projectionAlignmentRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.ProjectionAlignment;
		}
		if (sequenceAlignRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.Sequence;
		}
		if (fastaRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.Fasta;
		}
		if (loadMatrixRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.SqrMatrix;
		}
		if (numericalRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.Unknown;
		}
		if (geneticMapRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.GeneticMap;
		}
		if (tableReportRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.Table;
		}
		if (topmRadioButton.isSelected()) {
			return FileLoadPlugin.TasselFileType.TOPM;
		}
		return FileLoadPlugin.TasselFileType.Unknown;
	}

	public void okButton_actionPerformed(ActionEvent e) {
		isCancel = false;
		setVisible(false);
	}

	public void cancelButton_actionPerformed(ActionEvent e) {
		isCancel = true;
		setVisible(false);
	}

	public boolean isCancel() {
		return isCancel;
	}
}
