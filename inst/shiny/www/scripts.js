document.addEventListener("scroll", function(e) {
  const source = e.target;
  const bodies = "#tableOriginal .dataTables_scrollBody, #tableOriginal .dt-scroll-body, #tableMasked .dataTables_scrollBody, #tableMasked .dt-scroll-body";

  if(!(source instanceof Element) || !source.matches(bodies))
    return;

  const target = document.querySelector(source.closest("#tableOriginal") ?
    "#tableMasked .dataTables_scrollBody, #tableMasked .dt-scroll-body" :
    "#tableOriginal .dataTables_scrollBody, #tableOriginal .dt-scroll-body");

  if(target && target.scrollTop !== source.scrollTop)
    target.scrollTop = source.scrollTop;
}, true);

function markerDblclick(table) {
  table.on("dblclick", "tbody tr", function() {
    Shiny.setInputValue(
      "markerDblclick",
      table.row(this).index() + 1,
      {priority: "event"}
    );
  });
}

function sortBlanksLast(data, type) {
  if(type === "sort" && data === null) return Number.MAX_VALUE;
  if(type === "sort" && data === "") return "\uffff";
  return data;
}
