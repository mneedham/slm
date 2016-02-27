package org.neo4j.slm;

import java.util.ArrayList;
import java.util.List;

public class Node
{
    int nodeId;
    List<Relationship> in;
    List<Relationship> out;
    private int cluster;

    public Node( int nodeId )
    {
        this.nodeId = nodeId;
        this.in = new ArrayList<>();
        this.out = new ArrayList<>();
    }

    public Node in( Node source, double weight )
    {
        in.add( new Relationship( source.nodeId, this.nodeId, weight ) );
        return this;
    }

    public Node out( Node destination, double weight )
    {
        out.add( new Relationship( this.nodeId, destination.nodeId, weight ) );
        return this;
    }

    public int degree()
    {
        return out.size() + in.size();
    }

    public double weight() {
        double weight = 0.0;
        for ( Relationship relationship : in )
        {
            weight += relationship.getWeight();
        }

        for ( Relationship relationship : out )
        {
            weight += relationship.getWeight();
        }
        return weight;
    }

    public void setCluster( int cluster )
    {
        this.cluster = cluster;
    }

    public int getCluster()
    {
        return cluster;
    }



    public List<Relationship> relationships()
    {
        List<Relationship> rels = new ArrayList<>(  );
        for ( Relationship relationship : in )
        {
            rels.add(relationship);
        }
        for ( Relationship relationship : out )
        {
            rels.add(relationship);
        }
        return rels;
    }

    public List<Relationship> getIn()
    {
        return in;
    }

    public List<Relationship> getOut()
    {
        return out;
    }
}
