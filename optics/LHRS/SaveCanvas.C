//Author : Jin
//Save canvas as png, eps, cxx, root


#ifndef SAVE_CANVAS
#define SAVE_CANVAS


void SaveCanvas(TCanvas * c, TString name = "", Bool_t bEPS=kTRUE)
{
	if (name.Length()==0) name = c->GetName();

	TString rootfilename;

	c->Print(name+".png","png");
	c->Print(rootfilename = name+".root");

	if (bEPS)
	{
		gStyle->SetPaperSize( c->GetWindowWidth()/72*2.54,c->GetWindowHeight()/72*2.54);
		c->Print(name+".eps");
	}


}

#endif
