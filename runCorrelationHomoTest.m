

function[] = runCorrelationHomoTest(tissue)

tissue 
1
correlationHomoTest(1, 4000, 1, 4000, tissue);
tissue
2
correlationHomoTest(1, 4000, 4001, 8000, tissue);
tissue
3
correlationHomoTest(1, 4000, 8001, 12000, tissue);
tissue
4
correlationHomoTest(1, 4000, 12001, 18494, tissue);

tissue
5
correlationHomoTest(4001, 8000, 4001, 8000, tissue);
tissue 
6
correlationHomoTest(4001, 8000, 8001, 12000, tissue);
tissue
7
correlationHomoTest(4001, 8000, 12001, 18494, tissue);

tissue
8
correlationHomoTest(8001, 12000, 8001, 12000, tissue);
tissue 
9
correlationHomoTest(8001, 12000, 12001, 18494, tissue);

tissue
10
correlationHomoTest( 12001, 18494, 12001, 18494, tissue);
end