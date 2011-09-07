/*
 * ThreadedPluginListener
 */
package net.maizegenetics.plugindef;

/**
 *
 * @author terry
 */
public class ThreadedPluginListener extends Thread {

    private final PluginListener myPluginListener;
    private final PluginEvent myEvent;

    public ThreadedPluginListener(PluginListener pluginListener, PluginEvent event) {
        myPluginListener = pluginListener;
        myEvent = event;
    }

    public void run() {
        myPluginListener.dataSetReturned(myEvent);
    }

    public PluginListener getPluginListener() {
        return myPluginListener;
    }
}
