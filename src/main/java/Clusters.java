public class Clusters
{

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
}
