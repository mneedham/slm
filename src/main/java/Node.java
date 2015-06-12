import java.util.ArrayList;
import java.util.List;

public class Node
{
    int nodeId;
    List<ReducedNetwork.Relationship> in;
    List<ReducedNetwork.Relationship> out;
    private int cluster;

    public Node( int nodeId )
    {
        this.nodeId = nodeId;
        this.in = new ArrayList<>();
        this.out = new ArrayList<>();
    }

    public Node in( Node source, double weight )
    {
        in.add( new ReducedNetwork.Relationship( source.nodeId, this.nodeId, weight ) );
        return this;
    }

    public Node out( Node destination, double weight )
    {
        out.add( new ReducedNetwork.Relationship( this.nodeId, destination.nodeId, weight ) );
        return this;
    }

    public int degree()
    {
        return out.size() + in.size();
    }

    public double weight() {
        double weight = 0.0;
        for ( ReducedNetwork.Relationship relationship : in )
        {
            weight += relationship.weight();
        }

        for ( ReducedNetwork.Relationship relationship : out )
        {
            weight += relationship.weight();
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
}
