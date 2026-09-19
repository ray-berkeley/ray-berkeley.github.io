(function () {
  const body = document.body;
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const canvas = document.createElement("canvas");
  const peaks = [
    { x: 0.07, y: 0.16, rx: 68, ry: 24, angle: -0.20, phase: 0.3, depth: 18, levels: 6 },
    { x: 0.27, y: 0.29, rx: 42, ry: 15, angle: 0.28, phase: 1.7, depth: 10, levels: 5 },
    { x: 0.51, y: 0.12, rx: 30, ry: 12, angle: -0.45, phase: 3.1, depth: 14, levels: 4 },
    { x: 0.78, y: 0.19, rx: 58, ry: 20, angle: 0.12, phase: 4.4, depth: 22, levels: 6 },
    { x: 0.94, y: 0.43, rx: 48, ry: 17, angle: -0.32, phase: 2.4, depth: 16, levels: 5, negative: true },
    { x: 0.13, y: 0.69, rx: 52, ry: 19, angle: 0.40, phase: 5.2, depth: 20, levels: 5 },
    { x: 0.37, y: 0.84, rx: 66, ry: 23, angle: -0.12, phase: 3.8, depth: 12, levels: 6 },
    { x: 0.66, y: 0.66, rx: 38, ry: 14, angle: 0.34, phase: 0.9, depth: 18, levels: 5 },
    { x: 0.87, y: 0.85, rx: 62, ry: 21, angle: -0.38, phase: 2.9, depth: 24, levels: 6, negative: true }
  ];

  canvas.className = "ambient-contours";
  canvas.setAttribute("aria-hidden", "true");
  body.prepend(canvas);

  const context = canvas.getContext("2d");
  if (!context) {
    throw new Error("Canvas 2D context is unavailable");
  }

  let width = 0;
  let height = 0;
  let pointerTargetX = 0.5;
  let pointerTargetY = 0.5;
  let pointerX = 0.5;
  let pointerY = 0.5;
  let pointerActive = false;
  let lastFrame = 0;
  let revealed = false;

  function resizeCanvas() {
    const pixelRatio = Math.min(window.devicePixelRatio || 1, 2);
    width = window.innerWidth;
    height = window.innerHeight;
    canvas.width = Math.round(width * pixelRatio);
    canvas.height = Math.round(height * pixelRatio);
    context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);

    if (revealed) {
      drawContours(window.performance.now());
    }
  }

  function drawPeak(peak, time, opacityScale = 1) {
    const driftX = Math.sin(time * 0.000025 + peak.phase) * 8;
    const driftY = Math.cos(time * 0.000021 + peak.phase) * 6;
    const parallaxX = (pointerX - 0.5) * peak.depth;
    const parallaxY = (pointerY - 0.5) * peak.depth;
    const breath = 1 + Math.sin(time * 0.00005 + peak.phase) * 0.055;
    const color = peak.negative ? "190, 90, 70" : "0, 89, 152";

    context.save();
    context.translate(peak.x * width + driftX + parallaxX, peak.y * height + driftY + parallaxY);
    context.rotate(peak.angle + Math.sin(time * 0.000018 + peak.phase) * 0.035);
    context.lineWidth = 0.75;

    for (let level = 0; level < peak.levels; level += 1) {
      const scale = 1 - level * 0.12;
      const alpha = (0.065 + level * 0.008) * opacityScale;
      context.strokeStyle = `rgba(${color}, ${alpha})`;
      context.beginPath();
      context.ellipse(
        0,
        0,
        peak.rx * scale * breath,
        peak.ry * scale * (2 - breath),
        0,
        0,
        Math.PI * 2
      );
      context.stroke();
    }

    context.restore();
  }

  function drawContours(time) {
    context.clearRect(0, 0, width, height);
    peaks.forEach((peak) => drawPeak(peak, time));

    if (pointerActive) {
      drawPeak({
        x: pointerX,
        y: pointerY,
        rx: 34,
        ry: 13,
        angle: -0.25,
        phase: 1.2,
        depth: 0,
        levels: 4
      }, time, 0.45);
    }
  }

  function animate(time) {
    window.requestAnimationFrame(animate);
    if (time - lastFrame < 50) return;

    lastFrame = time;
    pointerX += (pointerTargetX - pointerX) * 0.06;
    pointerY += (pointerTargetY - pointerY) * 0.06;
    drawContours(time);
  }

  function followPointer(event) {
    pointerTargetX = event.clientX / window.innerWidth;
    pointerTargetY = event.clientY / window.innerHeight;
    pointerActive = true;
  }

  resizeCanvas();
  window.addEventListener("resize", resizeCanvas);

  if (!reducedMotion) {
    window.addEventListener("pointermove", followPointer, { passive: true });
  }

  window.setTimeout(function revealContours() {
    revealed = true;
    body.classList.add("ambient-contours-visible");
    drawContours(window.performance.now());

    if (!reducedMotion) {
      window.requestAnimationFrame(animate);
    }
  }, 30000);
})();
