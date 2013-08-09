package JISTIC;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

public class CombineMatricesNoParam {
	public static void main(String[] args) {
		String folderName = ".";
		if (args.length > 0)
			folderName = args[0];
		try {
			PrintStream saveoutput = System.out;

			File folder = new File(folderName);
			String[] files = folder.list();
			String[] alts={"AMP","DEL","LOH"};
			String[] ftype={"genes","miRNA"};
			for (int j = 0; j < alts.length; j++) {
				for (int k = 0; k < ftype.length; k++) {
					ArrayList<String> input = new ArrayList<String>();
					for (int i = 0; i < files.length; i++) {
						if (files[i].matches(alts[j]+"\\."+ftype[k]+"\\.chr.\\.matrix")
								|| files[i].matches(alts[j]+"\\."+ftype[k]+"\\.chr..\\.matrix")) {
							input.add(folderName + "/" + files[i] + " ");
						}
					}
					if (!input.isEmpty()) {
						System.setOut(new PrintStream(new FileOutputStream(folderName
								+ "/"+alts[j]+"."+ftype[k]+".All.matrix")));
						combine_matrices.main(input.toArray(new String[input.size()]));
						for (String string : input) {
							File f = new File(string);
							f.delete();							
						}
					}
				}
			}		
			System.setOut(saveoutput);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
