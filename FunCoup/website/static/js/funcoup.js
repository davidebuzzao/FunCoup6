// Functions for autocomplete, placed here, as it is used on multiple pages.
// All other javascript is put in the html-file where it is used

function autocomplete(inp, arr, type, badgeId, inputId) {
	/*the autocomplete function takes two arguments,
	the text field element and an array of possible autocompleted values:*/
	var currentFocus;

	/*execute a function when someone writes in the text field:*/
	inp.addEventListener("input", function(e) {
		showOptions(this)	
	});
	/*execute a function when someone click in the text field, only if not species box*/
	if(type!='speciesSelect'){
		inp.addEventListener("click", function(e) {
		showOptions(this)	
		});
	}
	/*execute a function presses a key on the keyboard:*/
	inp.addEventListener("keydown", function(e) {
		var x = document.getElementById(this.id + "autocomplete-list");
		if (x) x = x.getElementsByTagName("div");
		if (e.keyCode == 40) {
		  /*If the arrow DOWN key is pressed,
		  increase the currentFocus variable:*/
		  currentFocus++;
		  /*and and make the current item more visible:*/
		  addActive(x);
		} else if (e.keyCode == 38) { //up
		  /*If the arrow UP key is pressed,
		  decrease the currentFocus variable:*/
		  currentFocus--;
		  /*and and make the current item more visible:*/
		  addActive(x);
		} else if (e.keyCode == 13) {
		  /*If the ENTER key is pressed, prevent the form from being submitted,*/
		  e.preventDefault();
		  if (currentFocus > -1) {
			/*and simulate a click on the "active" item:*/
			if (x) x[currentFocus].click();
		  }
		}
	});

	function showOptions(t){
		var a, b, i, val = t.value;
		/*close any already open lists of autocompleted values*/
		closeAllLists();
		if (!val) { return false;}
		currentFocus = -1;
		/*create a DIV element that will contain the items (values):*/
		a = document.createElement("DIV");
		a.setAttribute("id", t.id + "autocomplete-list");
		a.setAttribute("class", "autocomplete-items");
		/*append the DIV element as a child of the autocomplete container:*/
		a.style.maxHeight = "300px";  // Limit the height to show 3-4 items
        a.style.overflowY = "auto";  // Enable vertical scrolling
		t.parentNode.appendChild(a);
		/*for each item in the array...*/
		if (arr){
            let count = 0;
            for (i = 0; i < arr.length && count < 4; i++) {  // Limit to 4 items
                var optionName = arr[i].name || arr[i];
                if (optionName.toUpperCase().includes(val.toUpperCase())) {
                    b = document.createElement("DIV");
                    var indexOfVal = optionName.toUpperCase().indexOf(val.toUpperCase());
                    b.innerHTML = optionName.substr(0, indexOfVal);
                    b.innerHTML += "<strong>" + optionName.substr(indexOfVal, val.length) + "</strong>";
                    b.innerHTML += optionName.substr(indexOfVal + val.length);
                    if (arr[i].url) {
                        b.innerHTML += "<input type='hidden' value='" + optionName + "' url='" + arr[i].url + "'>";
                    } else {
                        b.innerHTML += "<input type='hidden' value='" + optionName + "'>";
                    }
				
				/*execute a function when someone clicks on the item value (DIV element):*/
				b.addEventListener("click", function(e) {
					/*insert the value for the autocomplete text field:*/
					inp.value = this.getElementsByTagName("input")[0].value;
					if(type=='url'){
					window.location=this.getElementsByTagName("input")[0].getAttribute('url');
					}else if (type=='speciesSelect'){
						$('#'+badgeId).html(this.getElementsByTagName("input")[0].value+"<button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge("+"\""+badgeId+"\",\""+inputId+"\""+")'></button>")
						$("#"+inputId).css("color", "black");
						selectableSpecies();
						selectableCategories();
						selectableTissues();
					}else if(type=='filterPathway'){
						selectedPathway=this.getElementsByTagName("input")[0].value
						pathwayid=selectedPathway.split("pathwayID:")[1].replace(")","")
						pathwayName = selectedPathway.split("(pathwayID:")[0]
						existingHtml=$('#'+badgeId).html()
						$('#'+badgeId).html(existingHtml+"<div class='badge bg-secondary filter_badge text-wrap' id='"+pathwayid+"'>"+pathwayName+"<button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge("+"\""+pathwayid+"\",\"PATHWAYFILTER\""+")'></button></div>")
						inp.value="";
						// arr.splice(arr.indexOf(selectedPathway), 1); // Leaving commented as I dont know a way to put it back again...
						if (selectedPathway=="" || selectedPathway=="None"){
							p="None"
						}else{
							p=selectedPathway.split("pathwayID:")[1].replace(")","")
						}
						filterNetworkByPathway(p,"")

						// Check if the selected pathway is not already in the selectedPathways array
						if (!selectedFilterPathways.includes(selectedPathway)) {
							// Add the selected pathway to the selectedFilterPathways array
							selectedFilterPathways.push(selectedPathway);
						}
						// Update allFilterPathwayName by removing the selectedPathway
						allFilterPathwayName = allFilterPathwayName.filter(function(pathway) {
							return pathway !== selectedPathway;
						});
						allFilterPathwayName.sort();
						// Update the autocomplete dropdown with the filtered options
						autocomplete(document.getElementById("filterPathway"), allFilterPathwayName, 'filterPathway','pathwayFilterBadge','filterPathway');
					}else if(type=='filterTissue'){
						selectedTissue=this.getElementsByTagName("input")[0].value
						tissueID=selectedTissue.split("tissueID:")[1].replace(")","")
						tissueName = selectedTissue.split("(tissueID:")[0]
						existingHtml=$('#'+badgeId).html()
						$('#'+badgeId).html(existingHtml+"<div class='badge bg-secondary filter_badge text-wrap' id='"+tissueID+"'>"+tissueName+"<button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge("+"\""+tissueID+"\",\"TISSUEFILTER\""+")'></button></div>")
						inp.value="";
						// arr.splice(arr.indexOf(selectedTissue), 1);  // Leaving commented as I dont know a way to put it back again...
						filterNetworkByTissue(tissueID,"")

						// Check if the selected pathway is not already in the selectedPathways array
						if (!selectedFilterTissues.includes(selectedTissue)) {
							// Add the selected pathway to the selectedFilterTissues array
							selectedFilterTissues.push(selectedTissue);
						}
						// Update allFilterPathwayName by removing the selectedPathway
						allFilterTissueName = allFilterTissueName.filter(function(tissue) {
							return tissue !== selectedTissue;
						});
						allFilterTissueName.sort();
						// Update the autocomplete dropdown with the filtered options
						autocomplete(document.getElementById("filterTissue"), allFilterTissueName, 'filterTissue','tissueFilterBadge','filterTissue');
					}else if(type=='colorPathway'){
						// Update color and recreate badges
						updateColorAndBadges('colorPathway');
						selectedPathway=this.getElementsByTagName("input")[0].value
						pathwayid=selectedPathway.split("pathwayID:")[1].replace(")","")
						pathwayName=selectedPathway.split("(pathwayID:")[0]
						existingHtml=$('#'+badgeId).html()
						// TODO --> put text into one div and float-left; then put color picker and closing tab with float-right
						// $('#'+badgeId).html(existingHtml+"<div class='badge bg-secondary filter_badge text-wrap' id='"+pathwayid+"'>"+pathwayName+" <input type='color' class='colorPicker colorPickerPathway' id='pathwayColor_"+pathwayid+"' name='colorPicker' value='"+pathwayColor_dict[pathwayid]+"' onChange='changeNodeColorByAttribute(\""+pathwayid+"\",\"PATHWAYCOLOR\")'> <button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge(\""+pathwayid+"\",\"PATHWAYCOLOR\")'> </button> </div>")
						$('#'+badgeId).html(existingHtml+"<div class='badge bg-secondary filter_badge text-wrap' id='"+pathwayid+"'>"+pathwayName+" <input type='color' class='colorPicker colorPickerPathway' id='pathwayColor_"+pathwayid+"' name='colorPicker' value='"+pathwayColor_dict[pathwayid]+"' onChange='changeNodeColorByAttribute(\""+pathwayid+"\",\"PATHWAYCOLOR\")'> <button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge(\""+pathwayid+"\",\"PATHWAYCOLOR\")'> </button> </div>");

						inp.value="";
						// arr.splice(arr.indexOf(selectedPathway), 1); // Leaving commented as I dont know a way to put it back again...
						if (selectedPathway=="" || selectedPathway=="None"){
							p="None"
						}else{
							p=selectedPathway.split("pathwayID:")[1].replace(")","")
						}
						// Check if the selected pathway is not already in the selectedPathways array
						if (!selectedColorPathways.includes(selectedPathway)) {
							// Add the selected pathway to the selectedColorPathways array
							selectedColorPathways.push(selectedPathway);
						}
						// Update allFilterPathwayName by removing the selectedPathway
						allColorPathwayName = allColorPathwayName.filter(function(pathway) {
							return pathway !== selectedPathway;
						});
						allColorPathwayName.sort();

						if (!pathwayColor_badges.includes(pathwayid)) {
							// Add the selected pathway to the pathwayColor_badges array
							pathwayColor_badges.push(pathwayid)
						  }
						// Update the autocomplete dropdown with the filtered options
						autocomplete(document.getElementById("colorPathway"), allColorPathwayName, 'colorPathway','pathwayColorBadge','colorPathway');
					
					}else if(type=='colorTissue'){
						// Update color and recreate badges
						updateColorAndBadges('colorTissue');
						selectedTissue=this.getElementsByTagName("input")[0].value
						tissueID=selectedTissue.split("tissueID:")[1].replace(")","")
						tissueName = selectedTissue.split("(tissueID:")[0]
						existingHtml=$('#'+badgeId).html()
						$('#'+badgeId).html(existingHtml+"<div class='badge bg-secondary filter_badge text-wrap' id='"+tissueID+"'>"+tissueName+"  <input type='color' class='colorPicker colorPickerTissue' id='tissueColor_"+tissueID+"' name='colorPicker' value='"+tissueColor_dict[tissueID]+"' onChange='changeNodeColorByAttribute(\""+tissueID+"\",\"TISSUECOLOR\")'> <button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge(\""+tissueID+"\",\"TISSUECOLOR\")'> </button></div>")

						// $('#'+badgeId).html(existingHtml+"<div class='badge bg-secondary filter_badge text-wrap' id='"+tissueID+"'>"+tissueName+" <input type='color' class='colorPicker' id='tissueNodeColor' name='colorPicker' value=defaultNodeColor onChange='changeGlobalNodeColor()'/> <button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge("+"\""+tissueID+"\",\"TISSUECOLOR\""+")'></button></div>")
						inp.value="";
						// arr.splice(arr.indexOf(selectedTissue), 1);  // Leaving commented as I dont know a way to put it back again...
						// Check if the selected pathway is not already in the selectedPathways array
						if (!selectedColorTissues.includes(selectedTissue)) {
							// Add the selected pathway to the selectedColorTissues array
							selectedColorTissues.push(selectedTissue);
						}
						// Update allFilterPathwayName by removing the selectedPathway
						allColorTissueName = allColorTissueName.filter(function(tissue) {
							return tissue !== selectedTissue;
						});
						allColorTissueName.sort();

						if (!tissueColor_badges.includes(tissueID)) {
							// Add the selected pathway to the selectedFilterPathways array
							tissueColor_badges.push(tissueID)
						  }

						// Update the autocomplete dropdown with the filtered options
						autocomplete(document.getElementById("colorTissue"), allColorTissueName, 'colorTissue','tissueColorBadge','colorTissue');

					}else{
						$("#"+inputId).css("color", "black");
						$('#'+badgeId).html(this.getElementsByTagName("input")[0].value+"<button type='button' class='btn-close btn-close-white' aria-label='Close' onclick='removeBadge("+"\""+badgeId+"\",\""+inputId+"\""+")'></button>")
					}
				});
				a.appendChild(b);
				}
		}
		
		}

	}

	function addActive(x) {
	  /*a function to classify an item as "active":*/
	  if (!x) return false;
	  /*start by removing the "active" class on all items:*/
	  removeActive(x);
	  if (currentFocus >= x.length) currentFocus = 0;
	  if (currentFocus < 0) currentFocus = (x.length - 1);
	  x[currentFocus].classList.add("autocomplete-active");
	}
	
	function removeActive(x) {
	  /*a function to remove the "active" class from all autocomplete items:*/
	  for (var i = 0; i < x.length; i++) {
		x[i].classList.remove("autocomplete-active");
	  }
	}

	function closeAllLists(elmnt) {
		var x = document.getElementsByClassName("autocomplete-items");
		for (var i = 0; i < x.length; i++) {
			if (elmnt != x[i] && elmnt != inp) {
				x[i].parentNode.removeChild(x[i]);
	  		}
		}
  	}

	/*execute a function when someone clicks in the document:*/
	document.addEventListener("click", function (e) {
		if (e.target.parentNode && e.target.parentNode.className && e.target.parentNode.className!="autocomplete")
		closeAllLists(e.target);
	});

	inp.addEventListener("blur", function(e) {
        setTimeout(() => closeAllLists(), 200);  // Close the list when input loses focus, with a small delay to allow click to register
    });

}
