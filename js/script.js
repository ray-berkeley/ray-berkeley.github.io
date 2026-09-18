(function () {
  const root = document.documentElement;
  const body = document.body;
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  let targetX = 50;
  let targetY = 50;
  let currentX = 50;
  let currentY = 50;
  let animationFrame = null;

  function updateGradientPosition() {
    const deltaX = targetX - currentX;
    const deltaY = targetY - currentY;

    currentX += deltaX * 0.045;
    currentY += deltaY * 0.045;
    root.style.setProperty("--ambient-x", `${currentX.toFixed(2)}%`);
    root.style.setProperty("--ambient-y", `${currentY.toFixed(2)}%`);

    if (Math.abs(deltaX) + Math.abs(deltaY) > 0.02) {
      animationFrame = window.requestAnimationFrame(updateGradientPosition);
    } else {
      animationFrame = null;
    }
  }

  function followPointer(event) {
    targetX = (event.clientX / window.innerWidth) * 100;
    targetY = (event.clientY / window.innerHeight) * 100;

    if (animationFrame === null) {
      animationFrame = window.requestAnimationFrame(updateGradientPosition);
    }
  }

  if (!reducedMotion) {
    window.addEventListener("pointermove", followPointer, { passive: true });
  }

  window.setTimeout(function () {
    body.classList.add("ambient-gradient-visible");
  }, 30000);
})();
