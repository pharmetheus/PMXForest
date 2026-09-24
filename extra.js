document.addEventListener("DOMContentLoaded", () => {
  const articlesMenu = [...document.querySelectorAll(".navbar .dropdown-toggle")]
    .find(link => link.textContent.trim() === "Articles");

  if (articlesMenu) {
    articlesMenu.textContent = "Vignettes";
  }
});
