package org.panda.resource.signednetwork;

import org.biopax.paxtools.pattern.miner.SIFInteraction;
import org.panda.utility.CollectionUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ConflictResolver
{
	/**
	 * Inhibitions
	 */
	Set<SIFInteraction> inh;

	/**
	 * Positive phospho
	 */
	Set<SIFInteraction> pp;

	/**
	 * Positive phospho
	 */
	Set<SIFInteraction> np;

	/**
	 * Positive phospho
	 */
	Set<SIFInteraction> pe;

	/**
	 * Positive phospho
	 */
	Set<SIFInteraction> ne;

	/**
	 * Resolved as false.
	 */
	Set<SIFInteraction> removed;

	public ConflictResolver(Set<SIFInteraction> pp, Set<SIFInteraction> np, Set<SIFInteraction> pe,
		Set<SIFInteraction> ne, Set<SIFInteraction> inh)
	{
		this.pp = pp;
		this.np = np;
		this.pe = pe;
		this.ne = ne;
		this.inh = inh;
	}

	public Set<SIFInteraction> getPP()
	{
		return pp;
	}

	public Set<SIFInteraction> getNP()
	{
		return np;
	}

	public Set<SIFInteraction> getPE()
	{
		return pe;
	}

	public Set<SIFInteraction> getNE()
	{
		return ne;
	}

	public Set<SIFInteraction> getRemoved()
	{
		return removed;
	}

	public void decideAndRemoveConflictingInference()
	{
		Map<String, Set<String>> geneToIn = getGeneToInhibitors();
		removed = getRelationsToRemove(pp, np, geneToIn);
		removed.addAll(getRelationsToRemove(pe, ne, geneToIn));

		pp.removeAll(removed);
		np.removeAll(removed);
		pe.removeAll(removed);
		ne.removeAll(removed);
	}

	private Set<SIFInteraction> getRelationsToRemove(Set<SIFInteraction> pos, Set<SIFInteraction> neg,
		Map<String, Set<String>> geneToInh)
	{
		Set<SIFInteraction> rem = new HashSet<>();

		Map<String, SIFInteraction> posMap = convertToMap(pos);
		Map<String, SIFInteraction> negMap = convertToMap(neg);

		Set<String> common = new HashSet<>(posMap.keySet());
		common.retainAll(negMap.keySet());

		System.out.println("double edge size = " + common.size());

		for (String key : common)
		{
			String[] t = key.split(" ");
			String source = t[0];
			String target = t[1];

			if (geneToInh.containsKey(source))
			{
				for (String inh : geneToInh.get(source))
				{
					String iKey = inh + " " + target;

					if (posMap.containsKey(iKey))
					{
						rem.add(posMap.get(key));
					}
					else if (negMap.containsKey(iKey))
					{
						rem.add(negMap.get(key));
					}
				}
			}
		}

		System.out.println("resolved size = " + rem.size());

		return rem;
	}

	private Map<String, SIFInteraction> convertToMap(Set<SIFInteraction> sifs)
	{
		return sifs.stream().collect(Collectors.toMap(sif -> sif.sourceID + " " + sif.targetID, sif -> sif));
	}

	private Map<String, Set<String>> getGeneToInhibitors()
	{
		Map<String, Set<String>> map = new HashMap<>();
		for (SIFInteraction sif : inh)
		{
			String source = sif.sourceID;
			String target = sif.targetID;
			if (!map.containsKey(target)) map.put(target, new HashSet<>());
			map.get(target).add(source);
		}
		return map;
	}

	// Section: For investigating remaining conflicting relations.

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/PC/";
		String file = dir + "SignedREACH-woTF2.sif";
		Map<String, Map<String, Set<String>>> map = new HashMap<>();
		Files.lines(Paths.get(file)).map(l -> l.split("\t")).forEach(t ->
		{
			if (!map.containsKey(t[1])) map.put(t[1], new HashMap<>());
			String key = t[0] + " " + t[2];
			if (!map.get(t[1]).containsKey(key)) map.get(t[1]).put(key, new HashSet<>(Arrays.asList(t[3].split(";"))));
		});

		write(map, SignedType.PHOSPHORYLATES.getTag(), SignedType.DEPHOSPHORYLATES.getTag(), "phospho",
			dir + "conflict-phospho.sif");
		write(map, SignedType.UPREGULATES_EXPRESSION.getTag(), SignedType.DOWNREGULATES_EXPRESSION.getTag(), "expression",
			dir + "conflict-expression.sif");
	}

	private static void write(Map<String, Map<String, Set<String>>> map, String posType, String negType, String type,
		String filename) throws IOException
	{
		Set<String> overlap = new HashSet<>(map.get(posType).keySet());
		overlap.retainAll(map.get(negType).keySet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		for (String key : overlap)
		{
			Set<String> meds = new HashSet<>(map.get(posType).get(key));
			meds.addAll(map.get(negType).get(key));
			String[] t = key.split(" ");
			writer.write(t[0] + "\t" + type + "\t" + t[1] + "\t" + CollectionUtil.merge(meds, ";") + "\n");
		}

		writer.close();
	}
}
