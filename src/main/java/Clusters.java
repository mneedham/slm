import java.util.HashMap;
import java.util.Map;

public class Clusters
{
    private Map<Integer,Network.Cluster> clusters = new HashMap<>();

    public Network.Cluster get( int clusterId )
    {
        return clusters.get( clusterId );
    }

    public void put( int clusterId, Network.Cluster cluster )
    {
        clusters.put( clusterId, cluster );
    }

    public static int calculateNumberOfClusters( int[] cluster, int numberOfNodes )
    {
        if ( cluster == null )
        {
            return 0;
        }
        else
        {
            int i = 0;
            for ( int j = 0; j < numberOfNodes; j++ )
            {
                if ( cluster[j] > i )
                {
                    i = cluster[j];
                }
            }
            return i + 1;
        }
    }

    public int findClusterId( int nodeId )
    {
        for ( Map.Entry<Integer,Network.Cluster> entry : clusters.entrySet() )
        {
            for ( Node node : entry.getValue().nodes() )
            {
                if ( node.nodeId == nodeId )
                {
                    return entry.getKey();
                }
            }
        }

        return -1;
    }

    public boolean inSameCluster( int neighborNodeId, int nodeId )
    {
        int neighborCluster = -1;
        int nodeCluster = -1;

        for ( Map.Entry<Integer,Network.Cluster> entry : clusters.entrySet() )
        {
            for ( Node node : entry.getValue().nodes() )
            {
                if ( node.nodeId == nodeId )
                {
                    nodeCluster = entry.getKey();
                }

                if ( node.nodeId == neighborNodeId )
                {
                    neighborCluster = entry.getKey();
                }
            }
        }

        return neighborCluster != -1 && nodeCluster != -1 && (nodeCluster == neighborCluster);
    }

    public int[] nodesPerCluster( int numberOfNodes )
    {
        int[] numberOfNodesPerCluster = new int[numberOfNodes];
        for ( Map.Entry<Integer,Network.Cluster> entry : clusters.entrySet() )
        {
            numberOfNodesPerCluster[entry.getKey()] = entry.getValue().numberOfNodes();
        }
        return numberOfNodesPerCluster;
    }

    public int size()
    {
        return clusters.size();
    }
}
