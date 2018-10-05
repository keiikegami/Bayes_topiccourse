// Comparison

	decl mb_normal = loadmat("normal.csv");
	for(i=0;i<np;i++){
	DrawDensity(i,mb[][i]',sprint("Coef. No.",i+1),1,0);  DrawAdjust(ADJ_COLOR,2,10); 	DrawTitle(i, "");												
	DrawDensity(i,mb_normal[][i]',"Normal",1,0);  DrawAdjust(ADJ_COLOR,3,10); 	DrawTitle(i, "");	
	DrawAdjust(ADJ_AREA_X,i,-3,3,0);									
	}
	SaveDrawWindow(sprint("beta_comp.eps"));
  	CloseDrawWindow();