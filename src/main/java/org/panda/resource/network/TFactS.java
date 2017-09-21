package org.panda.resource.network;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

/**
 * Serves the TFactS database, the signed version, human, and has publication ID.
 *
 * To update this resource, get http://www.tfacts.org/TFactS-new/TFactS-v2/tfacts/data/Catalogues.xls, and save the
 * second sheet as tsv file.
 *
 * @author Ozgun Babur
 */
public class TFactS extends FileServer
{
	private static TFactS instance;

	private DirectedGraph positive;
	private DirectedGraph negative;
	private DirectedGraph unsigned;

	public static TFactS get()
	{
		if (instance == null) instance = new TFactS();
		return instance;
	}

	public DirectedGraph getPositiveGraph()
	{
		return positive;
	}

	public DirectedGraph getNegativeGraph()
	{
		return negative;
	}

	public DirectedGraph getUnsignedGraph()
	{
		return unsigned;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"TFactS_sign_sensitive_version2.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		positive = new DirectedGraph("TFactS positive", SignedType.UPREGULATES_EXPRESSION.getTag());
		negative = new DirectedGraph("TFactS negative", SignedType.DOWNREGULATES_EXPRESSION.getTag());
		unsigned = new DirectedGraph("TFactS all", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));
		while (scanner.hasNextLine())
		{
			String line = scanner.nextLine();

			// manual debug
			if (line.startsWith("RB1\tTP53\tDOWN\t")) continue;

			String[] token =  line.split("\t");

			if (token.length < 5 || !token[4].toLowerCase().contains("homo sapiens") || !hasPubmedID(token[5])) continue;

			switch (token[2])
			{
				case "UP":
					positive.putRelation(token[0], token[1]);
					break;
				case "DOWN":
					negative.putRelation(token[0], token[1]);
					break;
				default:
					System.err.println("Irregular sign: " + token[2]);
					break;
			}

			unsigned.putRelation(token[0], token[1]);
		}

		return true;
	}

	private boolean hasPubmedID(String s)
	{
		for (String t : s.split(";"))
		{
			if (!t.isEmpty())
			{
				try
				{
					int i = Integer.valueOf(t);
					if (i > 0) return true;
				}
				catch (NumberFormatException e) {}
			}
		}
		return false;
	}
	
	public static void main(String[] args)
	{
		boolean contains = get().getNegativeGraph().getDownstream("AATF").contains("CTNNB1");
		System.out.println("contains = " + contains);
	}
}
