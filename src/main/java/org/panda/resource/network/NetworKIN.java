package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.SiteSpecificGraph;
import org.panda.utility.statistics.UniquePrinter;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Serves the NetworKIN database.
 *
 * @author Ozgun Babur
 */
public class NetworKIN extends FileServer
{
	static NetworKIN instance;

	static SiteSpecificGraph graph;

	static Map<String, String> exceptionMapping;

	public static synchronized NetworKIN get()
	{
		if (instance == null) instance = new NetworKIN();
		return instance;
	}

	public SiteSpecificGraph getGraph()
	{
		return graph;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"NetworKIN.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://netphorest.info/download/networkin_human_predictions_3.1.tsv.xz"};
	}

	@Override
	public boolean load() throws IOException
	{
		UniquePrinter up = new UniquePrinter();
		graph = new SiteSpecificGraph("NetworKIN", SignedType.PHOSPHORYLATES.getTag());

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));

		// skip header
		scanner.nextLine();

		while (scanner.hasNextLine())
		{
			String line = scanner.nextLine();

			String[] token = line.split("\t");

			String source = token[2];
			// todo: map source to HGNC and remove the below ugly code
			source = source.replaceAll("alpha", "A").replaceAll("beta", "B")
				.replaceAll("delta", "D").replaceAll("eta", "H")
				.replaceAll("epsilon", "E").replaceAll("iota", "I")
				.replaceAll("theta", "Q").replaceAll("zeta", "Z")
				.replaceAll("III", "3").replaceAll("II", "2").toUpperCase();
			source = HGNC.get().getSymbol(source);

			if (source == null && exceptionMapping.containsKey(token[2]))
			{
				source = exceptionMapping.get(token[2]);
			}

			if (source == null)
			{
//				up.print("exceptionMapping.put(\"", token[2] + "\", \"\");");
				up.print("Not recognized: ", token[2]);
				continue;
			}

			String target = token[9];

			String targetUP = HGNC.get().getUniProt(target);

			if (targetUP == null) continue;

			int loc = Integer.valueOf(token[1]);

			String aa = UniProtSequence.get().getAminoacidAt(targetUP, loc);

			if (aa == null || !(aa.equals("S") || aa.equals("T") || aa.equals("Y"))) continue;

			String site = aa + loc;

			if (!graph.hasRelation(source, target))
			{
				graph.putRelation(source, target, "", site);
			}
			else
			{
				graph.addSite(source, target, site);
			}
		}

		return true;
	}

	public static void main(String[] args)
	{
		SiteSpecificGraph pn = get().getGraph();
//		pn.cropToNeighborhood(Collections.singleton("FLT3"));
//		pn.write("/home/ozgun/Documents/Temp/FLT3-phos.sif");

		pn.printStats();
	}

	static
	{
		exceptionMapping = new HashMap<>();
		exceptionMapping.put("PKCzeta", "PRKCZ");
		exceptionMapping.put("PKCdelta", "PRKCD");
		exceptionMapping.put("PKCeta", "PRKCH");
		exceptionMapping.put("PKCgamma", "PRKCG");
		exceptionMapping.put("PKCtheta", "PRKCQ");
		exceptionMapping.put("CK2a2", "CSNK2A2");
		exceptionMapping.put("PKAalpha", "PRKACA");
		exceptionMapping.put("PKAbeta", "PRKACB");
		exceptionMapping.put("PKAgamma", "PRKACG");
		exceptionMapping.put("CK2alpha", "CSNK2A1");
		exceptionMapping.put("CaMKIIgamma", "CAMK2G");
		exceptionMapping.put("CK1delta", "CSNK1D");
		exceptionMapping.put("CK1epsilon", "CSNK1E");
		exceptionMapping.put("CK1gamma2", "CSNK1G2");
		exceptionMapping.put("CK1gamma3", "CSNK1G3");
		exceptionMapping.put("PKBalpha", "AKT1");
		exceptionMapping.put("PKBbeta", "AKT2");
		exceptionMapping.put("PKBgamma", "AKT3");
		exceptionMapping.put("PKG1cGKI", "PRKG1");
		exceptionMapping.put("DMPK1", "DMPK");
		exceptionMapping.put("PDHK1", "PDK1");
		exceptionMapping.put("CaMKIV", "CAMK4");
		exceptionMapping.put("PKCepsilon", "PRKCE");
		exceptionMapping.put("ACTR2B", "ACVR2B");
		exceptionMapping.put("AuroraA", "AURKA");
		exceptionMapping.put("PDHK3", "PDK3");
		exceptionMapping.put("PDHK4", "PDK4");
		exceptionMapping.put("PKG2cGKII", "PRKG2");
		exceptionMapping.put("p70S6K", "RPS6KB1");
		exceptionMapping.put("AuroraB", "AURKB");
		exceptionMapping.put("AuroraC", "AURKC");
//		exceptionMapping.put("CK1alpha2", "");
		exceptionMapping.put("CaMKIalpha", "CAMK1");
		exceptionMapping.put("CaMKIdelta", "CAMK1D");
		exceptionMapping.put("CaMKIgamma", "CAMK1G");
	}
}
