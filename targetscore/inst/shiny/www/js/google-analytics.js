(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
})(window,document,'script','//www.google-analytics.com/analytics.js','ga');

ga('create', 'UA-57486113-13', 'auto');
ga('send', 'pageview');

$(document).ready(function() {
  $('input, textarea, select').change(function() {
    var tmp = {"db": "a", "p": "e", "items": ["b", "c", "d"]};
    console.log("ignore : " + JSON.stringify(tmp));
    ga('send', 'event', 'ignore', JSON.stringify(tmp));
  });
});
