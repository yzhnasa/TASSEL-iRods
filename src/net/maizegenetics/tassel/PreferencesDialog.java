package net.maizegenetics.tassel;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;

import net.maizegenetics.prefs.TasselPrefs;

/**
 * @author terryc
 */
public class PreferencesDialog extends JDialog {
    
    private final static Font HEADING_FONT = new Font(null, Font.BOLD, 14);
    private JCheckBox myRetainRareAlleles = new JCheckBox("Retain Rare Alleles");
    
    public PreferencesDialog() {
        super((Frame) null, "Preferences...", true);
        try {
            init();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private void init() throws Exception {
        
        setTitle("Preferences...");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        
        Container contentPane = getContentPane();
        
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.LEFT_ALIGNMENT);
        
        result.add(Box.createRigidArea(new Dimension(1, 30)));
        
        result.add(getGenotypeStoringPanel());
        
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        
        result.add(getButtons());
        
        result.add(Box.createRigidArea(new Dimension(1, 20)));
        
        contentPane.add(result);
        
        pack();
        
        setResizable(false);
        
    }
    
    private JPanel getGenotypeStoringPanel() {
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        
        JLabel alignPrefsLabel = new JLabel("Alignment Preferences...");
        alignPrefsLabel.setFont(HEADING_FONT);
        result.add(alignPrefsLabel);
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        
        myRetainRareAlleles.setSelected(TasselPrefs.getAlignmentRetainRareAlleles());
        result.add(myRetainRareAlleles);
        
        return result;
    }
    
    private JPanel getLine(String label, JTextField ref) {
        
        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));
        
        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);
        
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
    
    private void okButton_actionPerformed(ActionEvent e) {
        TasselPrefs.putAlignmentRetainRareAlleles(myRetainRareAlleles.isSelected());
        setVisible(false);
    }
    
    private void cancelButton_actionPerformed(ActionEvent e) {
        myRetainRareAlleles.setSelected(TasselPrefs.getAlignmentRetainRareAlleles());
        setVisible(false);
    }
}
