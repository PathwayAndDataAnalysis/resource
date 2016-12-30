package org.panda.resource.signednetwork;

import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.RelationshipXref;
import org.biopax.paxtools.model.level3.XReferrable;
import org.biopax.paxtools.model.level3.Xref;
import org.biopax.paxtools.pattern.miner.IDFetcher;
import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.biopax.paxtools.pattern.miner.SIFSearcher;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.panda.utility.Kronometre;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Generates the signed and directed version Pathway Commons SIF. This graph is produced to use in causality analysis.
 *
 * @author Ozgun Babur
 */
public class Generator
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/PC/";
		SimpleIOHandler h = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = h.convertFromOWL(new FileInputStream(dir + "PathwayCommons.8.Detailed.BIOPAX.owl"));
//		Model model = h.convertFromOWL(new FileInputStream(dir + "temp.owl"));
		System.out.println("Model size = " + model.getObjects().size());
//		BlacklistGenerator gen = new BlacklistGenerator();
//		Blacklist blacklist = gen.generateBlacklist(model);
//		blacklist.write(new FileOutputStream("../biopax-pattern/blacklist.txt"));
		Blacklist blacklist = new Blacklist(dir + "blacklist.txt");

		generate(model, blacklist, dir + "SignedPC.sif");
//		generate(model, blacklist, dir + "temp.sif");
	}

	public static void generate(Model model, Blacklist blacklist, String outFile) throws IOException
	{
		Kronometre kron = new Kronometre();

		SIFSearcher searcher;

		// prepare phospho-graph

		IDFetcher idFetcher = ele ->
		{
			if (ele instanceof XReferrable)
			{
				for (Xref xref : ((XReferrable) ele).getXref())
				{
					if (xref instanceof RelationshipXref && xref.getDb() != null &&
						xref.getDb().equalsIgnoreCase("hgnc symbol") && xref.getId() != null)
					{
						return Collections.singleton(xref.getId());
					}
				}
			}
			return Collections.emptySet();
		};

		searcher = new SIFSearcher(idFetcher, new PP1(), new PP2(), new PP3(), new PP4());
		searcher.setBlacklist(blacklist);
		Set<SignedSIFInteraction> pp = (Set<SignedSIFInteraction>) (Set<?>) searcher.searchSIF(model);

		System.out.println("Positive phospho = " + pp.size());

		searcher = new SIFSearcher(idFetcher, new PN1(), new PN2(), new PN3(), new PN4());
		searcher.setBlacklist(blacklist);
		Set<SignedSIFInteraction> pn = (Set<SignedSIFInteraction>) (Set<?>) searcher.searchSIF(model);

		System.out.println("Negative phospho = " + pn.size());

		decideConflictingPhosphorylation(pp, pn);

		// prepare expression graph

		searcher = new SIFSearcher(idFetcher, new EP1(), new EP2());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> ep = searcher.searchSIF(model);

		System.out.println("Positive expression = " + ep.size());

		searcher = new SIFSearcher(idFetcher, new EN1(), new EN2());
		searcher.setBlacklist(blacklist);
		Set<SIFInteraction> en = searcher.searchSIF(model);

		System.out.println("Negative expression " + en.size());

		decideConflictingExpression(ep, en);

		Set<SIFInteraction> sifs = new HashSet<SIFInteraction>();
		sifs.addAll(pp);
		sifs.addAll(pn);
		sifs.addAll(ep);
		sifs.addAll(en);

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (SIFInteraction sif : sifs)
		{
			writer.write(sif.toString(true) + "\n");
		}

		writer.close();
		kron.stop();
		kron.print();
	}

	private static void decideConflictingPhosphorylation(
		Set<SignedSIFInteraction> pos, Set<SignedSIFInteraction> neg)
	{
		Set<String> ov = getWOTypes(pos);
		ov.retainAll(getWOTypes(neg));
		System.out.println("\noverlap size = " + ov.size());


		Map<String, SignedSIFInteraction> posSIFs = mapSIFs(pos, ov);
		Map<String, SignedSIFInteraction> negSIFs = mapSIFs(neg, ov);

		assert posSIFs.size() == negSIFs.size();

		Set<String> posRem = new HashSet<String>();
		Set<String> negRem = new HashSet<String>();
		for (String key : posSIFs.keySet())
		{
			SignedSIFInteraction posSIF = posSIFs.get(key);
			SignedSIFInteraction negSIF = negSIFs.get(key);

			if (posSIF.p2med.isEmpty() && negSIF.mediators.size() > posSIF.mediators.size())
			{
				posRem.add(key);
				continue;
			}
			if (negSIF.p2med.isEmpty() && negSIF.mediators.size() < posSIF.mediators.size())
			{
				negRem.add(key);
				continue;
			}

			Set<String> phs = new HashSet<String>(posSIF.p2med.keySet());
			phs.retainAll(negSIF.p2med.keySet());

			Map<String, SignedSIFInteraction> decision = new HashMap<String, SignedSIFInteraction>();

			for (String ph : phs)
			{
				if (posSIF.p2med.get(ph).size() > negSIF.p2med.get(ph).size())
				{
					decision.put(ph, posSIF);
				}
				else if (posSIF.p2med.get(ph).size() < negSIF.p2med.get(ph).size())
				{
					decision.put(ph, negSIF);
				}
			}

			for (String ph : decision.keySet())
			{
				if (decision.get(ph) == posSIF)
				{
					negSIF.mediators.removeAll(negSIF.p2med.get(ph));
					negSIF.changedPhospho.remove(ph);
					negSIF.p2med.remove(ph);
				}
				else
				{
					posSIF.mediators.removeAll(posSIF.p2med.get(ph));
					posSIF.changedPhospho.remove(ph);
					posSIF.p2med.remove(ph);
				}
			}
		}

		removeDecided(pos, posRem);
		removeDecided(neg, negRem);

		System.out.println("decided = " + (posRem.size() + negRem.size()));
	}

	private static Map<String, SignedSIFInteraction> mapSIFs(Set<SignedSIFInteraction> sifs,
		Set<String> keys)
	{
		Map<String, SignedSIFInteraction> map = new HashMap<String, SignedSIFInteraction>();
		for (SignedSIFInteraction sif : sifs)
		{
			String key = sif.sourceID + "\t" + sif.targetID;
			if (keys.contains(key))
			{
				map.put(key, sif);
			}
		}
		return map;
	}

	private static void decideConflictingExpression(Set<SIFInteraction> pos, Set<SIFInteraction> neg)
	{
		Set<String> ov = getWOTypes(pos);
		ov.retainAll(getWOTypes(neg));
		System.out.println("\noverlap size = " + ov.size());


		Map<String, Integer> posScore = countMediators(pos, ov);
		Map<String, Integer> negScore = countMediators(neg, ov);

		assert posScore.size() == negScore.size();

		Set<String> posSet = new HashSet<String>();
		Set<String> negSet = new HashSet<String>();
		for (String key : posScore.keySet())
		{
			if (posScore.get(key) > negScore.get(key)) posSet.add(key);
			else if (posScore.get(key) < negScore.get(key)) negSet.add(key);
		}

		removeDecided(pos, negSet);
		removeDecided(neg, posSet);

		System.out.println("decided = " + (posSet.size() + negSet.size()));
	}

	private static Map<String, Integer> countMediators(Set<SIFInteraction> sifs, Set<String> st)
	{
		Map<String, Integer> map = new HashMap<String, Integer>();

		for (SIFInteraction sif : sifs)
		{
			String key = sif.sourceID + "\t" + sif.targetID;
			if (st.contains(key))
			{
				map.put(key, sif.mediators.size());
			}
		}
		return map;
	}

	private static void removeDecided(Set<? extends SIFInteraction> sifs, Set<String> keys)
	{
		Iterator<? extends SIFInteraction> iter = sifs.iterator();
		while (iter.hasNext())
		{
			SIFInteraction sif =  iter.next();
			String key = sif.sourceID + "\t" + sif.targetID;
			if (keys.contains(key)) iter.remove();
		}
	}

	private static String woType(SIFInteraction sif)
	{
		return sif.sourceID + "\t" + sif.targetID;
	}

	private static Set<String> getWOTypes(Set<? extends SIFInteraction> sifs)
	{
		Set<String> set = new HashSet<>(sifs.size());

		for (SIFInteraction sif : sifs)
		{
			set.add(woType(sif));
		}
		return set;
	}
}
