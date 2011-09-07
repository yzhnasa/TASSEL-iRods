package net.maizegenetics.baseplugins.alignment;

/**
 * Listener for detecting changes to a sequence.
 *
 * @author $Author: tcasstevens $
 * @version $Revision: 1.1 $, $Date: 2007/08/07 21:13:05 $ */
public interface SequenceListener
{
    public void sequenceNameChanged(Sequence aaseq, String old_name) throws Exception;
    public void sequenceAAChanged(Sequence aaseq);
    public void sequenceAnnotationChanged(Sequence aaseq);
    public void sequenceGroupChanged(Sequence aaseq);
    public void sequenceLineAnnotationsChanged(Sequence aaseq);
    public void sequenceColumnAnnotationsChanged(Sequence aaseq, int column);
    public void sequenceColorChanged(Sequence aasq);
}
