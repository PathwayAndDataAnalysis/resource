package org.panda.resource.signednetwork;

import org.panda.resource.SignedInteractionText;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Whenever double edge exists, this removes the negative ones.
 *
 * @author Ozgun Babur
 */
public class ConflictResolverSimple
{
	public void decideAndRemoveConflictingInference(String inFile, String resolvedFile, String removedFile) throws IOException
	{
		Set<SignedInteractionText> ints = Files.lines(Paths.get(inFile)).filter(l -> !l.isEmpty())
			.map(SignedInteractionText::new).collect(Collectors.toSet());

		// construct key to interaction mappings for different types

		Map<String, SignedInteractionText> ppMap = ints.stream().filter(i -> i.getType() == SignedType.PHOSPHORYLATES)
			.collect(Collectors.toMap(SignedInteractionText::signlessKey, i -> i));

		Map<String, SignedInteractionText> pnMap = ints.stream().filter(i -> i.getType() == SignedType.DEPHOSPHORYLATES)
			.collect(Collectors.toMap(SignedInteractionText::signlessKey, i -> i));

		Map<String, SignedInteractionText> epMap = ints.stream().filter(i -> i.getType() == SignedType.UPREGULATES_EXPRESSION)
			.collect(Collectors.toMap(SignedInteractionText::signlessKey, i -> i));

		Map<String, SignedInteractionText> enMap = ints.stream().filter(i -> i.getType() == SignedType.DOWNREGULATES_EXPRESSION)
			.collect(Collectors.toMap(SignedInteractionText::signlessKey, i -> i));

		Map<String, SignedInteractionText> gpMap = ints.stream().filter(i -> i.getType() == SignedType.ACTIVATES_GTPASE)
			.collect(Collectors.toMap(SignedInteractionText::signlessKey, i -> i));

		Map<String, SignedInteractionText> gnMap = ints.stream().filter(i -> i.getType() == SignedType.INHIBITS_GTPASE)
			.collect(Collectors.toMap(SignedInteractionText::signlessKey, i -> i));

		// map paired negative-positive relations to each other

		Map<SignedInteractionText, SignedInteractionText> pMap = pnMap.keySet().stream().filter(ppMap::containsKey)
			.collect(Collectors.toMap(pnMap::get, ppMap::get));

		Map<SignedInteractionText, SignedInteractionText> eMap = enMap.keySet().stream().filter(epMap::containsKey)
			.collect(Collectors.toMap(enMap::get, epMap::get));

		Map<SignedInteractionText, SignedInteractionText> gMap = gnMap.keySet().stream().filter(gpMap::containsKey)
			.collect(Collectors.toMap(gnMap::get, gpMap::get));

		Set<SignedInteractionText> remove = new HashSet<>();
		processAndSeparate(pMap, remove);
		processAndSeparate(eMap, remove);
		processAndSeparate(gMap, remove);

		ints.removeAll(remove);

		BufferedWriter w1 = Files.newBufferedWriter(Paths.get(resolvedFile));
		ints.forEach(i -> FileUtil.writeln(i.toString(), w1));
		w1.close();

		BufferedWriter w2 = Files.newBufferedWriter(Paths.get(removedFile));
		remove.forEach(i -> FileUtil.writeln(i.toString(), w2));
		w2.close();
	}

	void processAndSeparate(Map<SignedInteractionText, SignedInteractionText> map, Set<SignedInteractionText> remove)
	{
		map.keySet().stream().forEach(in ->
		{
			SignedInteractionText ip = map.get(in);
			in.subtractSites(ip);

			if (in.getSites() == null || in.getSites().isEmpty())
			{
				remove.add(in);
				in.addMediators(ip);
				in.setSites(in.getOrigSites());
			}
		});
	}

	public static void main(String[] args) throws IOException
	{
		ConflictResolverSimple crs = new ConflictResolverSimple();
		String dir = "/home/babur/Documents/PC/";
		crs.decideAndRemoveConflictingInference(dir + "SignedPC-dirty.sif", dir + "SignedPC.sif",
			dir + "SignedPC-removed.sif");
	}
}
