// Disable AJAX caching.
$.ajaxSetup({ cache: false });

function fetch_download_counts(){
  $.ajax({  url: "https://www.tbi.univie.ac.at/RNA/vrna_stats.json",
            type: "POST",
            dataType: "json",
            crossDomain: true,
            success: function(data, textStatus, jqXHR){
                      update_download_counters(data);
                     },
            error: function (responseData, textStatus, errorThrown) {
              console.warn(responseData, textStatus, errorThrown);
              alert('CORS failed - ' + textStatus);
            }
          });
}

function update_download_counters(data){
  for (var key in data) {
    if (data.hasOwnProperty(key)) {
      cid = ".dl-" + key.replace(/_/g, "-") + " span:first-child"

      $(cid).append("<span title=\"File Size\" class=\"sd-sphinx-override sd-badge sd-bg-warning sd-bg-text-warning\">" + data[key]['filesize']+ "</span>");

      if(data[key]['downloads'] >= 12)
        $(cid).append("<span title=\"Downloads\" class=\"sd-sphinx-override sd-badge sd-outline-muted sd-text-muted sd-rounded-pill float-end\">" + data[key]['downloads'] + "</span>");
    }
  }
}

$(document).ready(function(){
  setTimeout(fetch_download_counts, 2000);
});
