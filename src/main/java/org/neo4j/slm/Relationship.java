package org.neo4j.slm;

import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;

public class Relationship
{
    private final int source;
    private final int destination;
    private double weight;

    public Relationship( int source, int destination, double weight )
    {
        this.source = source;
        this.destination = destination;
        this.weight = weight;
    }

    static Relationship from( String line )
    {
        String[] splittedLine = line.split( "\t" );
        double weight = (splittedLine.length > 2) ? parseDouble( splittedLine[2] ) : 1;

        return new Relationship( parseInt( splittedLine[0] ), parseInt( splittedLine[1] ), weight );
    }

    public double getWeight()
    {
        return weight;
    }

    public int otherNode( int nodeId )
    {
        if ( nodeId == source )
        {
            return destination;
        }
        return source;
    }


    public int getSource()
    {
        return source;
    }

    public int getDestination()
    {
        return destination;
    }
}
