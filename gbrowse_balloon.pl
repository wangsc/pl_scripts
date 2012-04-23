 
balloon hover = sub   {
      my $feature = shift;
      my $box_content = "";
      if ($feature->strand > 0){
            $box_content .= "<b>Zoom in:</b>&nbsp;<a class=ggbmenu_link href=?q=NC_010473:"
               .($feature->start-25)."..".($feature->start+24)."&enable=DNA target=_self>".$feature->name." N-terminus</a><br>";
     }else{
            $box_content .= "<b>Zoom in:</b>&nbsp;<a class=ggbmenu_link href=?q=NC_010473:"
               .($feature->end-25)."..".($feature->end+24)."&enable=DNA target=_self>".$feature->name." N-terminus</a><br>";
     }
      return $box_content;
 }
