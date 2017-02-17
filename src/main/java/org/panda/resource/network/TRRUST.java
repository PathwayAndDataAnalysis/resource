package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.Graph;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

/**
 * Serves the TRRUST database.
 *
 * @author Ozgun Babur
 */
public class TRRUST extends FileServer
{
	private static TRRUST instance;

	private Graph unsigned;
	private Graph positive;
	private Graph negative;

	public static TRRUST get()
	{
		if (instance == null) instance = new TRRUST();
		return instance;
	}

	public Graph getUnsignedGraph()
	{
		return unsigned;
	}

	public Graph getPositiveGraph()
	{
		return positive;
	}

	public Graph getNegativeGraph()
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
		return new String[]{"http://www.grnpedia.org/trrust/trrust_rawdata.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		unsigned = new Graph("TRRUST all", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());
		positive = new Graph("TRRUST positive", SignedType.UPREGULATES_EXPRESSION.getTag());
		negative = new Graph("TRRUST negative", SignedType.DOWNREGULATES_EXPRESSION.getTag());

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));
		while (scanner.hasNextLine())
		{
			String[] token = scanner.nextLine().split("\t");
			if (token[2].equals("Activation")) positive.putRelation(token[0], token[1], true);
			if (token[2].equals("Repression")) negative.putRelation(token[0], token[1], true);
			unsigned.putRelation(token[0], token[1], true);
		}

		// Remove manually detected errors
		negative.removeRelation("ATM", "CDKN1A", true);

		return true;
	}

	public static void main(String[] args)
	{
		Graph upPC = SignedPC.get().getGraph(SignedType.UPREGULATES_EXPRESSION);
		Graph dwPC = SignedPC.get().getGraph(SignedType.DOWNREGULATES_EXPRESSION);
		Graph allPC = PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		Graph pos = get().getPositiveGraph();
		Graph neg = get().getNegativeGraph();
		Graph uns = get().getUnsignedGraph();

		allPC.printVennIntersections(uns);
		System.out.println();
		upPC.printVennIntersections(pos);
		System.out.println();
		dwPC.printVennIntersections(neg);
	}
}
