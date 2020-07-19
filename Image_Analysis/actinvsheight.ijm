//get input form the user
dir1 = getDirectory("Choose source Directory ");
list1 = getFileList(dir1);
n= list1.length;
height=newArray(n);
actin_top_area=newArray(n);
actin_bottom_area=newArray(n);
actin_top_tfi=newArray(n);
actin_bottom_tfi=newArray(n);
actin_total_tfi=newArray(n);
volume=newArray(n);
label=newArray(n);


for (i=0; i<list1.length; i++) {

		path = dir1+list1[i];
		open(path);
		run("8-bit");
		run("Clear Results");	
		//run("Threshold...");
		setAutoThreshold("IsoData dark");
		waitForUser("Dialog","Threshold for me");
		run("Set Measurements...", "area mean integrated stack limit display redirect=None decimal=4");
		setSlice(1);
		nslice=0;
		area=newArray(nSlices);
		tfi=newArray(nSlices);
		height[i]=0;
		for (j=0; j<nSlices; j++){
			setSlice(j+1);
			run("Measure");
			area[j]=getResult("Area",j);
			tfi[j]=getResult("RawIntDen",j);
			if (getResult("Area",j)>0){
				nslice=parseInt(getResult("Slice",j));
				height[i]=height[i]+0.5;
				}
			}
		
	
		low=parseInt(nslice-(height[i]/0.5));
		mid=parseInt(nslice-((height[i]/0.5)/2));
		actin_bottom_area[i]=0;
		actin_bottom_tfi[i]=0;
		actin_top_tfi[i]=0;
		actin_top_area[i]=0;
		actin_total_tfi[i]=0;
		volume[i]=0;

		for(j=parseInt(low); j<nslice; j++){
			if(j<mid){
				actin_bottom_area[i]=actin_bottom_area[i]+area[j];
				actin_bottom_tfi[i]=actin_bottom_tfi[i]+tfi[j];
				
				
			}
			else{
				actin_top_tfi[i]=actin_top_tfi[i]+tfi[j-1];
				actin_top_area[i]=actin_top_area[i]+area[j-1];
			}
			
		}
	
		actin_total_tfi[i]=actin_top_tfi[i] + actin_bottom_tfi[i];
		volume[i]=actin_top_area[i]+actin_bottom_area[i];

		label[i]=File.getName(path);
		
		
		close();

}

run("Clear Results");
for (k=0; k<list1.length; k++) {
	setResult("Label",k ,label[k]);
	setResult("Volume",k, volume[k]*0.5);
	setResult("Height",k,height[k]);
	setResult("Top_actin_Int",k,actin_top_tfi[k]);
	setResult("Bottom_actin_Int",k,actin_bottom_tfi[k]);
	setResult("Total_tfi",k,actin_total_tfi[k]);
	setResult("tfi_top/total", k, actin_top_tfi[k]/actin_total_tfi[k]);
	setResult("tfi_bottom/total",k,actin_bottom_tfi[k]/actin_total_tfi[k]);
	setResult("tfi_top/bottom",k,actin_top_tfi[k]/actin_bottom_tfi[k]);
	
	setResult("Top_actin_area",k,actin_top_area[k]);
	setResult("Bottom_actin_area",k,actin_bottom_area[k]);
	setResult("area_top/total", k, actin_top_area[k]/volume[k]);
	setResult("area_bottom/total",k,actin_bottom_area[k]/volume[k]);
	setResult("area_top/bottom",k,actin_top_area[k]/actin_bottom_area[k]);

	setResult("Total_actin_density",k,actin_total_tfi[k]/volume[k]);
	setResult("Top_actin_density",k,actin_top_tfi[k]/actin_top_tfi[k]);
	setResult("Bottom_actin_density",k,actin_bottom_tfi[k]/actin_bottom_tfi[k]);
	setResult("Top/Bottom density",k,(actin_top_tfi[k]*actin_bottom_area[k])/(actin_top_area[k]*actin_bottom_tfi[k]));
}
