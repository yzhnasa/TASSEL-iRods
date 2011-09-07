package net.maizegenetics.baseplugins.alignment;


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.HeadlessException;
import java.awt.Insets;
import java.util.Iterator;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.prefs.TasselPrefs;


/**
 * User: dkroon
 * Date: Mar 31, 2005
 * Time: 4:56:29 PM
 */
public class PreferencesDialog extends JDialog {
    
    private JPanel          pnlQualityScoreColors = new JPanel();
    private JPanel          pnlButtons = new JPanel();
    private JButton         btnOK;
    private JButton         btnCancel;
    private JButton         btnLowRangeColor;
    private JButton         btnMidRangeColor;
    private JButton         btnHighRangeColor;
    
    private JLabel          lblMidRange;
    
    private JSlider         lowRangeSlider;
    private JSlider         highRangeSlider;
    private JTextField      txtLowRange;
    private JTextField      txtMidRange;
    private JTextField      txtHighRange;
    private JPanel          pnlLowRangeColor;
    private JPanel          pnlMidRangeColor;
    private JPanel          pnlHighRangeColor;
    
    private Color           lowRangeColor;
    private Color           midRangeColor;
    private Color           highRangeColor;
    
    private int             lowRangeCutoff;
    private int             highRangeCutoff;
    private JTabbedPane     tabbedPane;
    private String          chooseColorButtonString         = "Set Color...";
    private String          jColorChooserString             = "Choose Quality Score Color";
    private JDialog         thisDialog;
    
    // This is not very good design.  Just trying to preserve
    // function until this package is refactored as a plugin.  -terry
    private static List <SequenceComponent> mySequenceComponents = new ArrayList();
    
    
    /**
     * Creates a non-modal dialog with the specified title and
     * with the specified owner frame.  If <code>owner</code>
     * is <code>null</code>, a shared, hidden frame will be set as the
     * owner of the dialog.
     * <p/>
     * This constructor sets the component's locale property to the value
     * returned by <code>JComponent.getDefaultLocale</code>.
     *
     * @param owner the <code>Frame</code> from which the dialog is displayed
     * @param title the <code>String</code> to display in the dialog's
     *              title bar
     * @throws java.awt.HeadlessException if GraphicsEnvironment.isHeadless()
     *                                    returns true.
     * @see java.awt.GraphicsEnvironment#isHeadless
     * @see javax.swing.JComponent#getDefaultLocale
     */
    public PreferencesDialog(Frame owner, String title) throws HeadlessException {
        super(owner, title);
        
        init();
        
        initGUI();
        
        initListeners();
    }
    
    
    private void init(){
        thisDialog      = this;
        
        int lrc = TasselPrefs.getAlignPluginQualscoreLowrangeColor();
        lowRangeColor = new Color(lrc);
        
        int mrc = TasselPrefs.getAlignPluginQualscoreMidrangeColor();
        midRangeColor = new Color(mrc);
        
        int hrc = TasselPrefs.getAlignPluginQualscoreHighrangeColor();
        highRangeColor  = new Color(hrc);
        
        int lowRangeCutoff = TasselPrefs.getAlignPluginQualscoreLowrangeCutoff();
        
        int highRangeCutoff = TasselPrefs.getAlignPluginQualscoreHighrangeCutoff();
    }
    
    private void initGUI() {
        this.getContentPane().setLayout(new BorderLayout());
        tabbedPane = new JTabbedPane();
        this.getContentPane().add(tabbedPane, BorderLayout.CENTER);
        pnlQualityScoreColors.setLayout(new GridBagLayout());
        
        btnOK = new JButton("OK");
        btnCancel = new JButton("Cancel");
        pnlButtons.add(btnOK);
        pnlButtons.add(btnCancel);
        this.getContentPane().add(pnlButtons, BorderLayout.SOUTH);
        
        lblMidRange = new JLabel("Middle Range Cutoff:");
        lblMidRange.setLabelFor(txtMidRange);
        pnlQualityScoreColors.add(lblMidRange, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        
        lowRangeSlider = new JSlider(-1, 99, lowRangeCutoff);
        lowRangeSlider.setBorder(BorderFactory.createTitledBorder("Low Range Cutoff:"));
        lowRangeSlider.setMajorTickSpacing(10);
        lowRangeSlider.setMinorTickSpacing(5);
        lowRangeSlider.setPaintTicks(true);
        lowRangeSlider.setPaintLabels(true);
        pnlQualityScoreColors.add(lowRangeSlider, new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 50, 0));
        
        highRangeSlider = new JSlider(-1, 99, highRangeCutoff);
        highRangeSlider.setBorder(BorderFactory.createTitledBorder("High Range Cutoff:"));
        highRangeSlider.setMajorTickSpacing(10);
        highRangeSlider.setMinorTickSpacing(5);
        highRangeSlider.setPaintTicks(true);
        highRangeSlider.setPaintLabels(true);
        pnlQualityScoreColors.add(highRangeSlider, new GridBagConstraints(0, 2, 3, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 50, 0));
        
        txtLowRange = new JTextField("" + lowRangeSlider.getValue());
        txtLowRange.setEditable(false);
        pnlQualityScoreColors.add(txtLowRange, new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 50, 0));
        
        txtMidRange = new JTextField("" + Math.min(lowRangeCutoff, highRangeCutoff) + " - " + Math.max(lowRangeCutoff, highRangeCutoff));
        txtMidRange.setEditable(false);
        pnlQualityScoreColors.add(txtMidRange, new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 50, 0));
        
        txtHighRange = new JTextField("" + highRangeSlider.getValue());
        txtHighRange.setEditable(false);
        pnlQualityScoreColors.add(txtHighRange, new GridBagConstraints(4, 2, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 50, 0));
        
        btnLowRangeColor = new JButton(chooseColorButtonString);
        pnlQualityScoreColors.add(btnLowRangeColor, new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));
        
        btnMidRangeColor = new JButton(chooseColorButtonString);
        pnlQualityScoreColors.add(btnMidRangeColor, new GridBagConstraints(6, 1, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));
        
        btnHighRangeColor = new JButton(chooseColorButtonString);
        pnlQualityScoreColors.add(btnHighRangeColor, new GridBagConstraints(6, 2, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));
        
        pnlLowRangeColor = new JPanel();
        pnlLowRangeColor.setBackground(lowRangeColor);
        pnlQualityScoreColors.add(pnlLowRangeColor, new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 10, 10));
        
        pnlMidRangeColor = new JPanel();
        pnlMidRangeColor.setBackground(midRangeColor);
        pnlQualityScoreColors.add(pnlMidRangeColor, new GridBagConstraints(7, 1, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 10, 10));
        
        pnlHighRangeColor = new JPanel();
        pnlHighRangeColor.setBackground(highRangeColor);
        pnlQualityScoreColors.add(pnlHighRangeColor, new GridBagConstraints(7, 2, 1, 1, 0.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 10, 10));
        
        tabbedPane.addTab("Quality Score Colors", pnlQualityScoreColors);
        tabbedPane.setSelectedComponent(pnlQualityScoreColors);
    }
    
    private void initListeners(){
        
        btnLowRangeColor.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                Color tempColor = JColorChooser.showDialog(thisDialog, jColorChooserString, lowRangeColor);
                if(tempColor != null){
                    lowRangeColor = tempColor;
                    pnlLowRangeColor.setBackground(lowRangeColor);
                }
            }
        });
        
        btnMidRangeColor.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                Color tempColor = JColorChooser.showDialog(thisDialog, jColorChooserString, midRangeColor);
                if(tempColor != null){
                    midRangeColor = tempColor;
                    pnlMidRangeColor.setBackground(midRangeColor);
                }
            }
        });
        
        btnHighRangeColor.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                Color tempColor = JColorChooser.showDialog( thisDialog, jColorChooserString, highRangeColor);
                if(tempColor != null){
                    highRangeColor = tempColor;
                    pnlHighRangeColor.setBackground(highRangeColor);
                }
            }
        });
        
        lowRangeSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                int value = lowRangeSlider.getValue();
                txtLowRange.setText("" + value);
                if(value > highRangeSlider.getValue()){
                    highRangeSlider.setValue(value);
                }
                setMidRangeTextField();
            }
        });
        
        
        
        highRangeSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                int value = highRangeSlider.getValue();
                txtHighRange.setText("" + value);
                if(value < lowRangeSlider.getValue()){
                    lowRangeSlider.setValue(value);
                }
                setMidRangeTextField();
            }
        });
        
        btnOK.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                TasselPrefs.putAlignPluginQualscoreLowrangeColor(lowRangeColor.getRGB());
                TasselPrefs.putAlignPluginQualscoreMidrangeColor(midRangeColor.getRGB());
                TasselPrefs.putAlignPluginQualscoreHighrangeColor(highRangeColor.getRGB());
                TasselPrefs.putAlignPluginQualscoreLowrangeCutoff(lowRangeSlider.getValue());
                TasselPrefs.putAlignPluginQualscoreHighrangeCutoff(highRangeSlider.getValue());
                
                Iterator itr = mySequenceComponents.iterator();
                while (itr.hasNext()) {
                    SequenceComponent current = (SequenceComponent) itr.next();
                    current.update();
                }
                
                thisDialog.setVisible(false);
            }
        });
        
        btnCancel.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                thisDialog.setVisible(false);
            }
        });
    }
    
    
    
    private void setMidRangeTextField(){
        int lowRange = lowRangeSlider.getValue();
        int highRange = highRangeSlider.getValue();
        if(lowRange == highRange ){
            txtMidRange.setText("N/A");
        }else{
            txtMidRange.setText("" + Math.min(lowRange, highRange) + " - " + Math.max(lowRange, highRange));
        }
    }
    
    public Color getLowRangeColor() {
        return lowRangeColor;
    }
    
    public Color getMidRangeColor() {
        return midRangeColor;
    }
    
    public Color getHighRangeColor() {
        return highRangeColor;
    }
    
    public int getLowRangeCutoff() {
        return lowRangeCutoff;
    }
    
    public int getHighRangeCutoff() {
        return highRangeCutoff;
    }
    
    /**
     * Returns the preferred size of this container.
     *
     * @return an instance of <code>Dimension</code> that represents
     *         the preferred size of this container.
     * @see #getMinimumSize
     * @see #getLayout
     * @see java.awt.LayoutManager#preferredLayoutSize(java.awt.Container)
     * @see java.awt.Component#getPreferredSize
     */
    public Dimension getPreferredSize() {
        return new Dimension(500,500);
    }
    
    
    public static void addSequenceComponent(SequenceComponent sc) {
        mySequenceComponents.add(sc);
    }
}
