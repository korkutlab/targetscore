Shiny.addCustomMessageHandler("showLoading",
	function(message) {
		console.log(message);
		
		if(message.show) {
			$("#loadingModal").modal("show");
		} else {
			$("#loadingModal").modal("hide");
		}
	}
);

Shiny.addCustomMessageHandler("resetFileInputs", function(message) {
	message.ids.forEach(function(id) {
		var fileInput = $("#" + id);
		var progress = $("#" + id + "_progress");

		fileInput.val("");
		fileInput.closest(".input-group").find('input[type="text"]').val("");
		progress.removeClass("active").css("visibility", "hidden");
		progress.find(".progress-bar")
			.removeClass("progress-bar-danger")
			.css("width", "0%")
			.text("");
		Shiny.setInputValue(id, null, {priority: "event"});
	});
});
