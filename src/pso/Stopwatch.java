/*=====================================================================
 *                      PlanWare Java Utilities
 *              Copyright (c) 2005, 2008 Kenneth B. Kopelson
 * 
 * Author  : Ken Kopelson
 * Created : Jun 30, 2005
 * File    : Stopwatch.java
 *
 *=====================================================================*/
package pso;

/**
 * A class for measuring the execution time of a program.
 * Allows the time to be tracked either in nanoseconds or in
 * milliseconds, depending on the setting useNanoTime constant.
 */
public class Stopwatch {
    public  static boolean useNanoTime  = false;
    private long    startTime    = -1L;
    private long    stopTime     = -1L;
    private long    restartTime  = -1L;
    private long    totalPaused  = 0L;
    private long    intervalTime = -1L;
    private boolean running   = false;

    /**
     * Starts the stopwatch, updating start time.
     * 
     * @return this stopwatch instance
     */
    public Stopwatch start() {
        startTime   = intervalTime = useNanoTime ? System.nanoTime() : System.currentTimeMillis();
        stopTime    = -1L;
        restartTime = -1L;
        totalPaused  = 0L;
        running     = true;
        return this;
    }
    
    /**
     * Stops the stopwatch, updating stop time.
     * 
     * @return this stopwatch instance
     */
    public Stopwatch stop() {
        if ( running ) { 
            stopTime = useNanoTime ? System.nanoTime() : System.currentTimeMillis();
            running  = false;
        }
        return this;
    }
    
    /**
     * Restarts this stopwatch, updating the total period of time for which it
     * was paused. Note that this method assumes the stopwatch had previously
     * been stopped.
     * 
     * @return this stopwatch instance
     */
    public Stopwatch restart() {
        if ( stopTime != -1L ) {
            restartTime = intervalTime = useNanoTime ? System.nanoTime() : System.currentTimeMillis();
            totalPaused += restartTime - stopTime;
            stopTime = -1L;
            running  = true;
        }
        return this;
    }
    
    /**
     * Gets the total elapsed time for this stopwatch, excluding any periods in
     * which it was paused.
     * 
     * @return total elapsed time
     */
    public long getElapsedTime() {
        if ( startTime == -1L )
            return 0;
        if ( running )
            return ( useNanoTime ? System.nanoTime() : System.currentTimeMillis() ) - startTime - totalPaused;
        else
            return stopTime - startTime - totalPaused;
    }
    
    /**
     * Returns the last interval, i.e. the last length of time for which
     * the stopwatch was running (if currently running, it will return the
     * time from the last start until now; if currently stopped, it will
     * return the time from last start until the last stop; if never run or
     * reset, it will return zero).
     * 
     * @return last interval
     */
    public long getLastInterval() {
        if ( startTime == -1L )
            return 0;
        if ( running ) {
            long curTime  = useNanoTime ? System.nanoTime() : System.currentTimeMillis();
            long interval = curTime - intervalTime;
            intervalTime  = curTime;
            return interval;
        }
        else
            return stopTime - intervalTime;
    }
    
    /**
     * Resets all time values for this stopwatch.
     * 
     * @return this stopwatch instance
     */
    public Stopwatch reset() {
        startTime    = -1L;
        intervalTime = -1L;
        stopTime     = -1L;
        restartTime  = -1L;
        totalPaused   = 0L;
        if ( running )
            start();
        return this;
    }
    
    /**
     * Returns the start time.
     * 
     * @return start time
     */
    public long getStartTime()   { return startTime; }
    
    /**
     * Returns the stop time.
     * 
     * @return stop time
     */
    public long getStopTime()    { return stopTime; }
    
    /**
     * Returns the total paused time.
     * 
     * @return total paused time
     */
    public long getTotalPaused() { return totalPaused; }
    
    /**
     * Checks whether stopwatch is currently running.
     * 
     * @return true if it is currently running, false otherwise
     */
    public boolean isRunning()   { return running; }
    
    /**
     * Returns a string representation of the current state
     * of the stopwatch, including elapsed time, length of
     * last interval, and total paused time. A custom message
     * can also be provided as an argument.
     * 
     * @param msg - The custom message
     * @return string representation
     */
    public String toString( String msg ) {
        long          interval = getLastInterval();
        long          paused   = getTotalPaused();
        StringBuilder sb       = new StringBuilder();
        String unitStr = useNanoTime ? "ns" : "ms";
        sb.append( "Elapsed: " ).append( getElapsedTime() ).append( unitStr );
        sb.append( " ( Interval: " ).append( interval ).append( unitStr );
        if ( paused > 0L )
            sb.append( ", Paused: " ).append( getTotalPaused() ).append( unitStr );
        sb.append( " )" );
        if ( msg != null )
            sb.append( ": " ).append( msg );
        return sb.toString();
    }

    @Override
	/**
	 * {@inheritDoc}
	 */
    public String toString() {
        return toString( null );
    }
    
    /**
     * Helper method for printing a message along
     * with a string representation of the stopwatch.
     * 
     * @param msg - The message to be printed
     */
    public void printMessage( String msg ) {
        System.out.println( toString( msg ));
    }
    
    /**
     * Helper method for printing the status of
     * the stopwatch without any additional messages.
     */
    public void printMessage() {
        printMessage( null );
    }
}
