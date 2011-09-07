package net.maizegenetics.tassel;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;
import net.maizegenetics.prefs.TasselPrefs;

/**
 * terryc
 */
public class PreferencesDialog extends JDialog {

    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton myStrictButton = new JRadioButton("Strict");
    private JRadioButton myNonStrictButton = new JRadioButton("Non-Strict");

    public PreferencesDialog() {
        super((Frame) null, "Strict...", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        setTitle("Strict...");
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

        myButtonGroup.add(myStrictButton);
        myButtonGroup.add(myNonStrictButton);

        if (TasselPrefs.getIDJoinStrict()) {
            myStrictButton.setSelected(true);
        } else {
            myNonStrictButton.setSelected(true);
        }

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

        String msg1 = "Choose Strict or Non-Strict...";
        String msg2 = "If Strict, taxa names must match exactly.";
        String msg3 = "If Non-Strict, taxa names must match up to";
        String msg4 = "least descriptive taxa.";

        result.add(Box.createRigidArea(new Dimension(1, 20)));
        result.add(new JLabel(msg1, JLabel.CENTER));
        result.add(new JLabel(msg2, JLabel.CENTER));
        result.add(new JLabel(msg3, JLabel.CENTER));
        result.add(new JLabel(msg4, JLabel.CENTER));
        result.add(Box.createRigidArea(new Dimension(1, 20)));

        //jLabel1.setFont(new Font("Dialog", Font.BOLD, 16));
        return result;

    }

    private JPanel getOptionPanel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        result.setBorder(BorderFactory.createEtchedBorder());

        result.add(myNonStrictButton);
        result.add(myStrictButton);

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

    public boolean isStrict() {
        return myStrictButton.isSelected();
    }

    private void okButton_actionPerformed(ActionEvent e) {
        if (myStrictButton.isSelected()) {
            TasselPrefs.putIDJoinStrict(true);
        } else {
            TasselPrefs.putIDJoinStrict(false);
        }
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        if (TasselPrefs.getIDJoinStrict()) {
            myStrictButton.setSelected(true);
        } else {
            myNonStrictButton.setSelected(true);
        }
        setVisible(false);
    }
}
