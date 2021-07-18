package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.FileUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.SiteSpecificGraph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Serves the TRRUST database.
 *
 * @author Ozgun Babur
 */
public class TRRUST extends FileServer
{
	private static TRRUST instance;

	private DirectedGraph unsigned;
	private DirectedGraph positive;
	private DirectedGraph negative;

	public static TRRUST get()
	{
		if (instance == null) instance = new TRRUST();
		return instance;
	}

	public DirectedGraph getUnsignedGraph()
	{
		return unsigned;
	}

	public DirectedGraph getPositiveGraph()
	{
		return positive;
	}

	public DirectedGraph getNegativeGraph()
	{
		return negative;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"TRRUST.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"};
	}

	@Override
	public boolean load() throws IOException
	{
		unsigned = new DirectedGraph("TRRUST all", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());
		positive = new DirectedGraph("TRRUST positive", SignedType.UPREGULATES_EXPRESSION.getTag());
		negative = new DirectedGraph("TRRUST negative", SignedType.DOWNREGULATES_EXPRESSION.getTag());

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));
		while (scanner.hasNextLine())
		{
			String[] token = scanner.nextLine().split("\t");
			Set<String> pmids = Arrays.stream(token[3].split(";")).map(s -> "PMID:" + s).collect(Collectors.toSet());
			if (token[2].equals("Activation")) positive.putRelation(token[0], token[1], pmids);
			if (token[2].equals("Repression")) negative.putRelation(token[0], token[1], pmids);
			unsigned.putRelation(token[0], token[1], pmids);
		}

//		DirectedGraph intersectingGraph = negative.getIntersectingGraph(positive);
//		negative.removeSubgraph(intersectingGraph);
//		positive.removeSubgraph(intersectingGraph);
//		intersectingGraph.printStats();

		// Remove manually detected errors
		negative.removeRelation("ATM", "CDKN1A");
		positive.removeRelation("RB1", "CDKN1A");

		return true;
	}

	public static void main(String[] args)
	{
//		Graph upPC = SignedPC.get().getGraph(SignedType.UPREGULATES_EXPRESSION);
//		Graph dwPC = SignedPC.get().getGraph(SignedType.DOWNREGULATES_EXPRESSION);
//		Graph allPC = PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
//		Graph pos = get().getPositiveGraph();
//		Graph neg = get().getNegativeGraph();
//		Graph uns = get().getUnsignedGraph();
//
//		allPC.printVennIntersections(uns);
//		System.out.println();
//		upPC.printVennIntersections(pos);
//		System.out.println();
//		dwPC.printVennIntersections(neg);

//		boolean contains = get().getNegativeGraph().getDownstream("TP53").contains("BCL2");
//		System.out.println("contains = " + contains);

		BufferedWriter writer = FileUtil.newBufferedWriter("/Users/ozgun/Documents/Data/PathwayCommonsV12/TRRust.sif");
		get().getPositiveGraph().write(writer);
		get().getNegativeGraph().write(writer);
		FileUtil.closeWriter(writer);

	}
}
