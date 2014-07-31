/*
 * FileLoadPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.trait.ReadPhenotypeUtils;
import net.maizegenetics.dna.snp.ReadPolymorphismUtils;
import net.maizegenetics.dna.snp.ReadSequenceAlignmentUtils;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.util.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.prefs.TasselPrefs;

import org.apache.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import net.maizegenetics.dna.map.TOPMUtils;

/**
 *
 * @author Ed Buckler
 */
public class FileLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FileLoadPlugin.class);
    private String[] myOpenFiles = null;
    private TasselFileType myFileType = TasselFileType.Unknown;
    private PlinkLoadPlugin myPlinkLoadPlugin = null;
    private ProjectionLoadPlugin myProjectionLoadPlugin = null;
    private JFileChooser myOpenFileChooser = new JFileChooser(TasselPrefs.getOpenDir());

    public enum TasselFileType {

        SqrMatrix, Sequence, Unknown, Fasta,
        Hapmap, Plink, Phenotype, ProjectionAlignment, Phylip_Seq, Phylip_Inter, GeneticMap, Table,
        Serial, HapmapDiploid, Text, VCF, HDF5, TOPM
    };
    public static final String FILE_EXT_HAPMAP = ".hmp.txt";
    public static final String FILE_EXT_HAPMAP_GZ = ".hmp.txt.gz";
    public static final String FILE_EXT_PLINK_MAP = ".plk.map";
    public static final String FILE_EXT_PLINK_PED = ".plk.ped";
    public static final String FILE_EXT_SERIAL_GZ = ".serial.gz";
    public static final String FILE_EXT_HDF5 = ".h5";
    public static final String FILE_EXT_VCF = ".vcf";
    public static final String FILE_EXT_TOPM = ".topm";
    public static final String FILE_EXT_TOPM_H5 = ".topm.h5";
    public static final String FILE_EXT_FASTA = ".fasta";

    /**
     * Creates a new instance of FileLoadPlugin
     */
    public FileLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public FileLoadPlugin(Frame parentFrame, boolean isInteractive, PlinkLoadPlugin plinkLoadPlugin, ProjectionLoadPlugin projectionLoadPlugin) {
        super(parentFrame, isInteractive);
        myPlinkLoadPlugin = plinkLoadPlugin;
        myProjectionLoadPlugin = projectionLoadPlugin;
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
                
                setOpenFiles(getOpenFilesByChooser());
                theDialog.dispose();
            }

            if ((myOpenFiles == null) || (myOpenFiles.length == 0)) {
                return null;
            }

            List result = new ArrayList();
            ArrayList<String> alreadyLoaded = new ArrayList();
            for (int i = 0; i < myOpenFiles.length; i++) {

                if (alreadyLoaded.contains(myOpenFiles[i])) {
                    continue;
                }

                DataSet tds = null;
                try {

                    if (myFileType == TasselFileType.Unknown) {
                        if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP) || myOpenFiles[i].endsWith(FILE_EXT_HAPMAP_GZ)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Hapmap);
                        } else if ((myOpenFiles[i].endsWith(FILE_EXT_TOPM_H5)) || (myOpenFiles[i].endsWith(FILE_EXT_TOPM))) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.TOPM);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.TOPM);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_PED)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                            String theMapFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_PED, FILE_EXT_PLINK_MAP);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(theMapFile);
                            myPlinkLoadPlugin.loadFile(myOpenFiles[i], theMapFile, null);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_MAP)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                            String thePedFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_MAP, FILE_EXT_PLINK_PED);
                            alreadyLoaded.add(myOpenFiles[i]);
                            alreadyLoaded.add(thePedFile);
                            myPlinkLoadPlugin.loadFile(thePedFile, myOpenFiles[i], null);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_SERIAL_GZ)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Serial);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Serial);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_HDF5)) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.HDF5);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.HDF5);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_VCF) || myOpenFiles[i].endsWith(FILE_EXT_VCF + ".gz")) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.VCF);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.VCF);
                        } else if (myOpenFiles[i].endsWith(FILE_EXT_FASTA) || myOpenFiles[i].endsWith(FILE_EXT_FASTA + ".gz")) {
                            myLogger.info("guessAtUnknowns: type: " + TasselFileType.Fasta);
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = processDatum(myOpenFiles[i], TasselFileType.Fasta);
                        } else {
                            alreadyLoaded.add(myOpenFiles[i]);
                            tds = guessAtUnknowns(myOpenFiles[i]);
                        }
                    } else {
                        alreadyLoaded.add(myOpenFiles[i]);
                        tds = processDatum(myOpenFiles[i], myFileType);
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    StringBuilder builder = new StringBuilder();
                    builder.append("Error loading: ");
                    builder.append(myOpenFiles[i]);
                    builder.append("\n");
                    builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
                    String str = builder.toString();
                    if (isInteractive()) {
                        DialogUtils.showError(str, getParentFrame());
                    } else {
                        myLogger.error(str);
                    }
                }

                if (tds != null) {
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, FileLoadPlugin.class));
                }

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
                    while (guess == null && inline != null && (inline.startsWith("#") || inline.startsWith("<"))) {
                        if (inline.startsWith("<")) {
                            String[] info = tagPattern.split(inline);
                            if (info[1].toUpperCase().startsWith("MARKER")) {
                                isMarker = true;
                            } else if (info[1].toUpperCase().startsWith("TRAIT")) {
                                isTrait = true;
                            } else if (info[1].toUpperCase().startsWith("NUMER")) {
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
                    throw new IOException("Improperly formatted header. Data will not be imported.");
                }
            } else if ((line1.startsWith(">")) || (line1.startsWith(";"))) {
                guess = TasselFileType.Fasta;
            } else if (sval1.length == 1) {
                guess = TasselFileType.SqrMatrix;
            } else if ((line1.startsWith("#Nexus")) || (line1.startsWith("#NEXUS")) || (line1.startsWith("CLUSTAL"))
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

    private DataSet processDatum(String inFile, TasselFileType theFT) {
        Object result = null;
        String suffix = null;
        try {
            switch (theFT) {
                case Hapmap: {
                    suffix = FILE_EXT_HAPMAP;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_HAPMAP_GZ;
                    }
                    result = ImportUtils.readFromHapmap(inFile, this);
                    break;
                }
                case HDF5: {
                    IHDF5Reader reader=HDF5Factory.openForReading(inFile);
                    boolean t4HDF5=HDF5Utils.isTASSEL4HDF5Format(HDF5Factory.openForReading(inFile));
                    reader.close();
                    if(t4HDF5) {
                        DialogUtils.showError("This file is TASSEL4 HDF5 format. It will be converted to TASSEL5 " +
                                "HDF5 format with the t5.h5 suffix added.  This may take a few minutes.",getParentFrame());
                        String newInfile=inFile.replace(".h5",".t5.h5");
                        MigrateHDF5FromT4T5.copyGenotypes(inFile,newInfile);
                        inFile=newInfile;
                    }
                    suffix = FILE_EXT_HDF5;
                    result = ImportUtils.readGuessFormat(inFile);
                    break;
                }
                case VCF: {
                    suffix = FILE_EXT_VCF;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_VCF + ".gz";
                    }
                    result = ImportUtils.readFromVCF(inFile, this);
                    break;
                }
                case Sequence: {
                    result = ReadSequenceAlignmentUtils.readBasicAlignments(inFile, 40);
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
                    result = TableReportUtils.readDelimitedTableReport(inFile, "\t");
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
            builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
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
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        return null;
    }

    /**
     * Provides a open filer that remember the last location something was
     * opened from
     */
    private File[] getOpenFilesByChooser() {
        myOpenFileChooser.setMultiSelectionEnabled(true);
        File[] lopenFiles = null;
        myOpenFileChooser.setVisible(true);
        int returnVal = myOpenFileChooser.showOpenDialog(getParentFrame());
        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            lopenFiles = myOpenFileChooser.getSelectedFiles();
            String test = myOpenFileChooser.getCurrentDirectory().getPath();
            TasselPrefs.putOpenDir(myOpenFileChooser.getCurrentDirectory().getPath());
        }
        return lopenFiles;
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

    public void setOpenFiles(String[] openFiles) {
        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
        } else {
            myOpenFiles = openFiles;
        }
    }

    public TasselFileType getTheFileType() {
        return myFileType;
    }

    public void setTheFileType(TasselFileType theFileType) {
        myFileType = theFileType;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FileLoadPlugin.class.getResource("/net/maizegenetics/analysis/images/LoadFile.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Load";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load data from files on your computer.";
    }
}

/**
 * <p>Title: TASSEL</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: USDA-ARS</p>
 *
 * @author Edward Buckler
 * @version 1.0
 */
class FileLoadPluginDialog extends JDialog {

    boolean isCancel = true;
    ButtonGroup conversionButtonGroup = new ButtonGroup();
    JRadioButton hapMapRadioButton = new JRadioButton("Load Hapmap");
    JRadioButton hdf5RadioButton = new JRadioButton("Load HDF5");
    JRadioButton vcfRadioButton = new JRadioButton("Load VCF");
    JRadioButton plinkRadioButton = new JRadioButton("Load Plink");
    JRadioButton sequenceAlignRadioButton = new JRadioButton("Load Phylip");
    JRadioButton fastaRadioButton = new JRadioButton("Load FASTA File");
    JRadioButton numericalRadioButton = new JRadioButton("Load Numerical (trait, covariates, or factors)");
    JRadioButton loadMatrixRadioButton = new JRadioButton("Load Square Numerical Matrix (i.e. kinship)");
    JRadioButton guessRadioButton = new JRadioButton("Make Best Guess");
    JRadioButton projectionAlignmentRadioButton = new JRadioButton("Load Projection Alignment");
    JRadioButton geneticMapRadioButton = new JRadioButton("Load a Genetic Map");
    JRadioButton tableReportRadioButton = new JRadioButton("Load a Table Report");
    JRadioButton topmRadioButton = new JRadioButton("Load a TOPM (Tags on Physical Map)");

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
