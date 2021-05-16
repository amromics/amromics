<style scoped>
</style>
<template>
<div id="piechart">
</div>
</template>
<script>
/* eslint-disable */
import {PanPieChart} from "@/amromicsjs";
import SampleAPI from "@/api/SampleAPI";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'PangenomePieChart',
    props: ['core_data_url'],
    data() {
        return {
            loading: false,
            core_data:undefined
        };
    },
    computed: {
    collectionId() {
      return this.$route.params.cid;
      ;
    }
    },
    async mounted() {
      this.loading = true;
      var ctx=document.getElementById('piechart');
    //console.log(this.core_data);
      var pie = new PanPieChart(ctx);
      const value = await SampleAPI.fetchPangenomeSummary(this.collectionId);
      this.core_data=value.data;
      console.log(this.core_data.group);
     
      pie.load(this.core_data.group);
      
      pie.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
